all:
	markdown README.md > README.html
	doxygen > /dev/null 2> /dev/null
