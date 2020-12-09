(wget "http://www.genesilico.pl/modomics/api/modifications?id=1&format=csv" -q -O modomics.csv ) && 
	sed -i -e "s:</br>:\n:g" -e "s:<br>:\n:g" -e "s:<pre>::g" -e "s:</pre>::g" modomics.csv &&
	cat modomics.csv | python trim_lines.py >> modomics2.csv && 
	rm modomics.csv &&
	mv modomics2.csv modomics.csv
