freq_vector = [3 8 15 30 35 40] 

sine_matrix = [];

t = 0:0.001:3

i = 1; 

for freq_iterator = freq_vector

	temp = [];
	temp = sin(2*pi*freq_iterator*t);
	sine_matrix(:,i) = temp ;
	i = i+1;
end	
