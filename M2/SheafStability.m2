-*
Copyright 2019 Jonathan Steinbuch

You may redistribute this file under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3 of
the License, or any later version.
*-

newPackage(
           "SheafStability",
           Version=>"0.1.0",
           Date=>"August 19, 2019",
           Authors=>{{Name=> "Jonathan Steinbuch",
                    Email=>"jonathan.steinbuch@uni-osnabrueck.de"}},
           Headline=>"Macaulay2 package to use a function to compute sheaf stability",
           DebuggingMode => false,
	   Certification => {
		"journal name" => "",
		"journal URI" => "",
		"article title" => "",
		"acceptance date" => "",
		"published article URI" => "",
		"published code URI" => "",
		"repository code URI" => "",
		"release at publication" => 0,	    -- as an integer
		"version at publication" => "",
		"volume number" => "",
		"volume URI" => ""
		}
)

export{"computeSemistability"}

computeSemistability = method(TypicalValue => Boolean)
computeSemistability Matrix := M -> (
	S := ring M;
	I := gens gb ideal presentation S;
	R := ring I;
	
	infilename := "input.txt";
	outfilename := "output_test2.txt"; 
	
	infile := infilename << "";
	infile << "semistability" << endl;
	infile << "characteristic: " << 0 << endl;
	infile << "variables: ";
	Rvars := flatten entries vars R;
	for i in 0 .. (length Rvars)-2 do (
	infile << "\"" << toString (Rvars#i) << "\"" << ", ";
	);
	infile  << "\"" << toString (Rvars#((length Rvars)-1)) << "\"";
	infile << endl;
	rels := flatten entries I;
	infile << "relations: ";
		for i in 0 .. (length rels)-2 do (
	infile << toString (rels#i) << ", ";
	);
	infile << toString (rels#((length rels)-1));
	infile << endl;
	infile << "matrix: " << toString entries M << endl;
	infile << close;
	
	cmd := concatenate("../Debug/stability --input-file=\"", infilename, "\" --output-file=\"",outfilename,"\"");
 	run cmd;
 	
 	outfile := get outfilename;
 	outlist := lines outfile;
 	
 	retVal := true;
 	
 	if value(outlist#0) == 1 then (
 		Ker := value(outlist#2);
 		table := value(concatenate("new HashTable from ", outlist#1));
 		retVal = (false,table,Ker);
 	) else (
 		retVal = true;
 	);
 
 	
 	return retVal;
);


-*  *-