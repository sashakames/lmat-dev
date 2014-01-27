{ len=length($(base)); 

mink=30;  

if (len<50) {
    mink=len-20;   
}
show=1; 
if  ( ($NF == "DirectMatch" || $NF == "MultiMatch" ) && ($(base+3)>= mink))  { 

    for (i = base + 5 ; i < (NF-2) ; i+=2) { 

	if ($(i) == 9606) 
	{ 

	    show=0;
	    break;
	} 
    } 
    if (show == 1) 
	print ">" $1 " " $2 " " $(NF -2) " " $(base+3) "\n" $(base);
#	print;
} }

