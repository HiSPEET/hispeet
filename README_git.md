     mkdir HiSPEET
     cd HiSPEET

     git init
     mkdir doc
     mkdir src
     ..
     git add doc src
     git commit

     # stage the modified and deleted files
		 git add -u

		 # add and commit only the modified and deleted files.
		 git commit -a
