	- warum verwenden wir eher die elastic net modelle?
		- weil wir für die elastik net modelle eine performance estimation haben
		- und weil wir weniger features als samples haben (bei ols maximale anzahl an features ist die Anzahl ansamples (low rank problem )- 			  somit ist das ein limiteirender faktor - wobei wir bei mir sowiso wenig features haben)

	- warum is es unser Ziel möglichst viele features/ regionen pro Modell zu haben?
		- einfach höhere Chance das eine Region dabei ist, welche wirklich gut ist?
		- Laura: unser Ziel ist es viele regionene zufinden um diese weiter untersuchen zu können - mehr Signal = höhere Chance was zu finden
			- bei Ihec daten sind es beispielsweise 300 - 400 regionen pro gen
			- wenn wir viel finden, erwarte wir das wir auch wirklich etwas mit der expression  zu tun hat - viele features welche danne
			  influss auf die expression haben
	
	- es kann sein das der correlation estimate (the output coefficient? - nein) auf den 5 genen des referenzsets von eg 5 nicht brauchbar und somit 
	  die p valiue auch nicht brauchbar 
		--> damit sind die berechnete p-values basierend auf den der corss validation gemeint