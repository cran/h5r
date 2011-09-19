
export RCMD=R-devel

clean : 
	rm -f *.Rout *.h5 .Rhistory

all :
	${RCMD} CMD BATCH --no-save testRead.R &
	${RCMD} CMD BATCH --no-save testReadPoints.R &
	${RCMD} CMD BATCH --no-save testWrite.R &
	${RCMD} CMD BATCH --no-save testErrors.R & 
	${RCMD} CMD BATCH --no-save testGetGroup.R &
	${RCMD} CMD BATCH --no-save testCompound.R & 

valgrind:
	${RCMD} -d "valgrind --tool=memcheck --leak-check=full" --vanilla < testGetGroup.R


