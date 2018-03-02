require 'irb'
require 'pp'

NPROC = [1,2,4,8,16,32,64,128,256,512,1024]
SIZES = [2048, 3072, 4096, 6144, 8192, 12288, 16384]

def parse_file(fn)
	r=/--- Приведение к треугольному виду: (\S+) сек\n--- Обратный ход метода Гаусса: (\S+)/
	s=File.read fn
	if s =~ r
		{triag:$1.to_f, gauss:$2.to_f}
	else
		nil
	end
end

def parse_dir(dir)
	hsh = {}

	NPROC.each do |nproc|
		hsh[nproc] = {}
		SIZES.each do |size|
			hsh[nproc][size] = parse_file("#{dir}/#{size}_#{nproc}")
		end
	end
	# pp hsh
	hsh
end

def av(hsh1,hsh2)
	hsh3 = {}

	NPROC.each do |nproc|
		hsh3[nproc] = {}
		SIZES.each do |size|
			if !hsh1[nproc][size].nil?
				hsh3[nproc][size] = {}
				hsh3[nproc][size][:triag] = (hsh1[nproc][size][:triag]+hsh2[nproc][size][:triag])/2
				hsh3[nproc][size][:gauss] = (hsh1[nproc][size][:gauss]+hsh2[nproc][size][:gauss])/2
			else
				hsh3[nproc][size]=nil
			end
		end
	end
	hsh3
end

def print_all(hsh)
	NPROC.each do |nproc|
		SIZES.each do |size|
			if !hsh[nproc][size].nil?
				printf "#{hsh[nproc][size][:triag].round(6)}"
			else
				printf "Timeout!"
			end
			printf ", " if size != SIZES.max
		end
		puts
	end
	puts
	puts
	NPROC.each do |nproc|
		SIZES.each do |size|
			if !hsh[nproc][size].nil?
				printf "#{hsh[nproc][size][:gauss].round(6)}"
			else
				printf "Timeout!"
			end
			printf ", " if size != SIZES.max
		end
		puts
	end
end

def parse_dirs(dir_list)
	hsh = parse_dir dir_list[0]

	if dir_list.size >= 2
		dir_list[1..-1].each do |dir|
			hsh = av(hsh, parse_dir(dir))
		end
	end
	hsh
end

dir_list = []
base = 'full'
(1..2).each do |i|
	dir_list << base+i.to_s
end

hsh = parse_dirs dir_list
print_all hsh