.PHONY: doc

doc:
	rm -rf doc/html
	ford doc/hispeet.md

distclean:
	rm -rf build*
	rm -rf doc/html
	find . -name .DS_Store -exec rm {} \;
