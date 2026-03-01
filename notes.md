Without traditional matrix-vector multiplication, problem set 1 takes 30.73 
seconds

With MKL, it takes 40.51 seconds. But of course it's slower, since I
am having to allocate, deallocate, and copy several arrays each time


RefineOneVar uses compute_Az