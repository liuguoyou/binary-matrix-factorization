for lmi in 0 1 2 3; 
do 
	for mi in 0 1 2 3 4;
	do
		for du in 0 1; 
		do 
			src/bsvd_test -l 6 -k 800 -d $du -L $lmi -m 1 -M 1 -i ${mi}   data/mnist_data.pbm  > salida_du${du}_lmi${lmi}_mi${mi}.txt 
			mv atoms_mosaic.pbm atoms_mosaic_du${du}_lmi${lmi}_mi${mi}.pbm
			mv residual_mosaic.pbm residual_mosaic_${du}_lmi${lmi}_mi${mi}.pbm
		done 
	done
done
