all: ROSBot0.wrl

ROSBot0.wrl: ROSBot.py EZCAD3.py
	ROSBot.py

doxygen:
	markdown README.md > README.html
	doxygen > /dev/null 2> /dev/null
