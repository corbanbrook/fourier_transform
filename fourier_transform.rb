require 'benchmark'
require 'optparse'

options = {}
options[:samplerate] = 44100.0
options[:waveform] = :sine
options[:benchmark] = false
options[:transform] = :fft

option_parser = OptionParser.new do |opts|
  opts.banner = "Usage: ruby fourier_transform.rb -n buffersize -f frequency [-r samplerate] [-w sine|triangle|square|saw] [--plot] [--benchmark]"
  opts.on("-n", "--buffersize N", Integer, "Specify buffer size of N samples") do |n|
    options[:buffersize] = n
  end
  opts.on("-f", "--frequency FREQUENCY", "Specify the frequency of the generated signal in Hz") do |f|
    options[:frequency] = f.to_f
  end
  opts.on("-r", "--samplerate SAMPLERATE", Float, "Specify the sample rate of the generated signal") do |r|
    options[:samplerate] = r
  end
  opts.on("-w", "--waveform WAVEFORM", "Specify the shape of the waveform. sine (default), triangle, square, or saw") do |wave|
    options[:waveform] = wave.to_sym if ['sine', 'triangle', 'square', 'saw'].include?(wave)
  end
  opts.on("--dft", "Perform a Discrete Fourier Transform (slower)") do
    options[:transform] = :dft
  end
  opts.on("--plot", "Plot a graph of the Fourier Transform") do |plot|
    options[:plot] = plot
  end
  opts.on("--benchmark", "Benchmark DFT vs FFT") do |b|
    options[:benchmark] = b
  end
  opts.on_tail("-h", "--help", "Show this message") do
      puts opts
      exit
  end
end

begin
  option_parser.parse!
rescue
  puts $! # print out error
  option_parser.parse('--help') # print out command glossary
end

class Float
  def round_to n = 0
    (self * 10**n).round / 10.0**n
  end
end

class FourierTransform
  attr_reader :spectrum, :bandwidth, :samplerate, :buffersize
  
  def initialize buffersize, samplerate
    @buffersize = buffersize
    @samplerate = samplerate
    @bandwidth = (2.0 / @buffersize) * (@samplerate / 2.0)
    @spectrum = Array.new
    
    build_reverse_table
    build_trig_tables
  end

  def build_reverse_table
    @reverse = Array.new(@buffersize)
    @reverse[0] = 0;

    limit = 1
    bit = @buffersize >> 1

    while (limit < @buffersize )
      (0...limit).each do |i|
        @reverse[i + limit] = @reverse[i] + bit
      end
    
      limit = limit << 1
      bit = bit >> 1
    end
  end

  def build_trig_tables
    @sin_lookup = Array.new(@buffersize)
    @cos_lookup = Array.new(@buffersize)
    (0...@buffersize).each do |i|
      @sin_lookup[i] = Math.sin(- Math::PI / i);
      @cos_lookup[i] = Math.cos(- Math::PI / i);
    end
  end
  
  def dft(buffer)
    real = Array.new(buffer.length/2, 0)
    imag = Array.new(buffer.length/2, 0)
    
    (0...buffer.length/2).each do |k|
      (0...buffer.length).each do |n|
        real[k] += buffer[n] * Math.cos(2 * Math::PI * k * n / buffer.length)
        imag[k] += buffer[n] * -Math.sin(2 * Math::PI * k * n / buffer.length)
      end
      @spectrum[k] = 2 * Math.sqrt(real[k] ** 2 + imag[k] ** 2) / buffer.length
    end
    
    @spectrum
  end

  def fft(buffer)
    raise Exception if buffer.length % 2 != 0 

    real = Array.new(buffer.length)
    imag = Array.new(buffer.length)

    (0...buffer.length).each do |i|
      real[i] = buffer[@reverse[i]]
      imag[i] = 0.0
    end

    # here begins teh Danielson-Lanczos section
    halfsize = 1
    while halfsize < buffer.length
      #k = - Math::PI / halfsize
      #phase_shift_step_real = Math.cos(k)
      #phase_shift_step_imag = Math.sin(k)
      phase_shift_step_real = @cos_lookup[halfsize]
      phase_shift_step_imag = @sin_lookup[halfsize]
      current_phase_shift_real = 1.0
      current_phase_shift_imag = 0.0
      (0...halfsize).each do |fft_step|
        i = fft_step
        while i < buffer.length
          off = i + halfsize
          tr = (current_phase_shift_real * real[off]) - (current_phase_shift_imag * imag[off])
          ti = (current_phase_shift_real * imag[off]) + (current_phase_shift_imag * real[off])
          real[off] = real[i] - tr
          imag[off] = imag[i] - ti
          real[i] += tr
          imag[i] += ti

          i += halfsize << 1
        end
        tmp_real = current_phase_shift_real
        current_phase_shift_real = (tmp_real * phase_shift_step_real) - (current_phase_shift_imag * phase_shift_step_imag)
        current_phase_shift_imag = (tmp_real * phase_shift_step_imag) + (current_phase_shift_imag * phase_shift_step_real)
      end

      halfsize = halfsize << 1
    end

    (0...buffer.length/2).each do |i|
      @spectrum[i] = 2 * Math.sqrt(real[i] ** 2 + imag[i] ** 2) / buffer.length
    end
    
    @spectrum
  end

  def index_to_frequency(i)
    i * @bandwidth
  end
  
  def frequency_to_index(freq)
    (@buffersize.to_f * (freq / @samplerate)).round
  end
  
  def peak_frequency
    index = (0...spectrum.length).max {|a, b| spectrum[a] <=> spectrum[b] }
    index_to_frequency(index)
  end
  
  def plot(rows = 20, cols = 80)
    return if @spectrum.empty?

    max = @spectrum.max
    min = @spectrum.min
    y = (max - min) / rows.to_f
    bandwidth = cols / @spectrum.length
    rows.downto(0).each do |row|
      line = ""
      (0...@spectrum.length).each do |col|
        if row == 0
          line << "-"
        elsif @spectrum[col].round_to(1) >= (row * y).round_to(1)
          line << "."
        else
          line << " "
        end
      end
      if row % 2 == 0
        line << "- #{sprintf "%.1f", row * y}"
      end

      puts line
    end

    # Draw Freq arrows
    line = ""
    (0..@spectrum.length).each do |i|
      if @spectrum[i] == max
        line << "*"
      elsif i % 10 == 0
        line << "^"
      else 
        line << " "
      end
    end
    puts line  

    # Draw Frequency labels
    line = ""
    (0..@spectrum.length).each do |col|
      if col % 10 == 0
        label = "#{(col * @bandwidth / 1000).round}kHz"
        if col == 0
          line << label
        else
          line = line.chop << label + " "
        end
        (0...(10-label.length)).each do 
          line << " "
        end
      end
    end
    puts line  
  end
end

# generate a signal to pass to the FourierTransform
signal = Array.new
(0...options[:buffersize]).each do |i|
  step = i * (options[:frequency] / options[:samplerate].to_f) % 1.0
  case options[:waveform]
    when :sine # no harmonics
      signal[i] = Math.sin(step * 2 * Math::PI) 
    when :square # lots of harmonics
      signal[i] = step < 0.5 ? 1.0 : -1.0
    when :triangle
      signal[i] = 1 - 4 * (step.round - step).abs
    when :saw
      signal[i] = 2 * (step - step.round)
  end
end

fourier = FourierTransform.new(options[:buffersize], options[:samplerate])

if options[:benchmark]
  Benchmark.bm(10) do |bench|
    bench.report('DFT') do
      fourier.dft(signal)
    end
    bench.report('FFT') do
      fourier.fft(signal)
    end
  end
else
  fourier.send(options[:transform], signal) # runs either fft or dft transform

  puts "[#{options[:transform].to_s.upcase}] Sample rate: #{fourier.samplerate/1000}kHz / Buffer size: #{fourier.buffersize} samples / Input: generated #{options[:frequency]}Hz #{options[:waveform].to_s} wave\n\n"
  puts "      Found fundamental peak frequency of #{fourier.peak_frequency.round_to(2)}Hz +/- #{(fourier.bandwidth/2.0).round_to(2)} (off by #{(options[:frequency] - fourier.peak_frequency).round_to(2).abs}Hz)\n\n"  

  if options[:plot]
    fourier.plot
  end
end