FLAGS  = -c # -CB -fpe0 -traceback#O5   
PROF =  
LIB =   libacml.a lib2decomp_fft.a
COMP  = mpif90 -132 -r8  -mcmodel=medium -shared_intel 

post:  decomp_2d.o post.o stuff.o par_post.txt
	$(COMP) decomp_2d.o io.o post.o stuff.o $(LIB) -o bud
post.o : post.f par_post.txt 
	$(COMP) $(FLAGS) post.f
#ders_n.o : ders_n.f par_post.txt
#	$(COMP) $(FLAGS) ders_n.f
stuff.o : stuff.f par_post.txt
	$(COMP) $(FLAGS) stuff.f
#spectral.o : spectral.f par_post.txt
#	$(COMP) $(FLAGS) spectral.f
#fftpack5.o : fftpack5.f90 par_post.txt
#	$(COMP) $(FLAGS) fftpack5.f90
decomp_2d.o: decomp_2d.f90
	$(COMP)  -DDOUBLE_PREC -DOVERWRITE -O3  -cpp  -c decomp_2d.f90 io.f90 
clean:
	rm *.o *.mod bud
