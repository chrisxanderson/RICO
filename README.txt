RICO: Random Input Correlated Output

There are two main parts to this project.
1) The rico R-module
2) The RicoDoc.pdf file (and supporting LaTeX objects)

## To create the rico module ##
Remain in this top-level directory execute (%> is a shell prompt)
%> make

This will create a file called "rico_0.1.tar.gz" (module the version number). Create a
shell script to load rico and run R. Here's what I use,

-- in a Bash Shell -------------------------------------------------------------
#!/usr/bin/env sh
export R_DEFAULT_PACKAGES='datasets,utils,grDevices,graphics,stats,methods,rico'
R --quiet --no-save --no-restore
--------------------------------------------------------------------------------

To test if rico is running in R type (at the R> prompt)
> x <- rico.norm()
> plot(x)

The result should be plot of a standard normal probability density function.

## Documentation ##
The RicoDoc.pdf file is located in ./Documentation/RicoDoc.pdf
The R-module is supposed to come with lots of documentation. 
Perhaps some is there as you read this.

## Other Files and Directories ##
./Library          : contains PDF's and saved  webpages of relevant external documentation
./notes            : contains some code snippets used to develop the R module and
./notes/notes.txt  : contains is the developers log book for creating the R module. 
		   : The notes are chronological top-to-bottom
                   : If all else fails, the answer may be in this file. Start from the bottom.
./TODO.txt         : The TODO items for the R-module. There is a separate TODO for RicoDoc

