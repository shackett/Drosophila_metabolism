import re
import sys

class ReactionCompound:
    """Represents a compound combined with its stoichiometry"""
    def __init__(self, name, number=1):
        """
        Initialized with the name of the compound and optionally the number
        of molecules it has in a formula
        """
        self.name = name
        self.number = number

    def __str__(self):
        if self.number > 1:
            return str(self.number) + " " + self.name
        else:
            return "\"" + self.name + "\""


class Enzyme:
    """Represents an enzyme and the reaction that it catalyzes"""
    def __init__(self, name, reversible=True):
        self.name = name
        self.reactants = []
        self.products = []
        self.reversible = reversible

    def add_reactant(self, reactant):
        self.reactants.append(reactant)

    def add_product(self, product):
        self.products.append(product)
    
    def __str__(self):
        ret = "\"" + self.name + "\": "
        ret += " + ".join(map(str, self.reactants))
        if self.reversible:
            ret += " <=> "
        else:
            ret += " -> "
        ret += " + ".join(map(str, self.products))
        return ret

mycompound = ReactionCompound("pseudochemicaline", 3)

myenzyme = Enzyme("pseudoenzymine")
myenzyme.add_reactant(ReactionCompound("A", 3))
myenzyme.add_reactant(ReactionCompound("B"))
myenzyme.add_product(ReactionCompound("C", 2))

myenzyme2 = Enzyme("pseudoenzymine")
myenzyme2.add_reactant(ReactionCompound("A", 3))
myenzyme2.add_reactant(ReactionCompound("B"))
myenzyme2.add_product(ReactionCompound("C", 4))
######## write a system of reaction equations to be read by Omix using reaction output from R ######

Rreactions = open('meta.rxns.tsv','r')

enzymes = []
for line in Rreactions:
	rxn = re.compile('reaction')
	react = re.compile('reactant')
	prod = re.compile('product')
	
	terms = line.split("\t")
	if rxn.search(line):
		if terms[2] == '1':
			reversible = False
		else:
			reversible = True
		enzymes.append(Enzyme(terms[1], reversible))
		
	elif react.search(line):
		current_enzyme = enzymes.pop()
		current_enzyme.add_reactant(ReactionCompound("".join(["\"", terms[1], "\"", "[", terms[3][0], "]"]), terms[2]))
		enzymes.append(current_enzyme)
		
	elif prod.search(line):
		
		current_enzyme = enzymes.pop()
		current_enzyme.add_product(ReactionCompound("".join(["\"", terms[1], "\"", "[", terms[3][0], "]"]), terms[2]))
		enzymes.append(current_enzyme)
		
rxn_eqtn = open('reaction_equations','w')		
for enzyme in enzymes:
	rxn_eqtn.write("%s\n" % (enzyme))
rxn_eqtn.close()

	
	
	
	
	