function [final] = Calculate_W_matrix(w, dist)
    global ph_size grid_size r n_transducers c wn f_start f_end fh_start fh_end wnh snr;


    W_matrix = dist;
    final = zeros(1,length(dist));
    
    for i=1:length(dist)
       d = W_matrix(i);
       x = (w/c)*d;
       z = complex(0, x);
       temp = (w * exp(z)) / d;
       j = complex(0, -1);
       temp = temp * j;
       final(i) = temp;
    end
    final = reshape(final, 1, []);
end

