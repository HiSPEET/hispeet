.PHONY: doc

doc:
	rm -rf doc/html
	ford doc/hispeet.md

distclean:
	rm -rf build*
	rm -rf bin/__pycache__
	rm -rf doc/html
