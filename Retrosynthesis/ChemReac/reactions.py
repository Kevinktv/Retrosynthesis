"""
This file contains various reactions
"""
import re

from rdkit import Chem
from rdkit.Chem import AllChem, Draw

from Retrosynthesis.ChemReac.functional_groups import FunctionalGroups


############################
def draw(smile):
    """
    Given a smile, this function draws it.
    :param smile: The smile to draw
    """
    molecule = Chem.MolFromSmiles(smile)

    fig = Draw.MolToMPL(molecule, size=(200, 200))


############################
class Reactions:

    @staticmethod
    def isAromatic(smile):
        # Checks whether the smile is aromatic or not. I am assuming that a compound is aromatic if u got a cycle of benzene.
        # So a molecule is aromatic if it contains a ring of 6 carbons. This is denoted
        # in smiles as two small letter c1.
        return smile.count('c1') == 2 and smile.count('c') == 6

    @staticmethod
    def isAlkene(smile):
        # Return whether smile has a C=C bond.
        return False if smile.find('C=C') == -1 else True

    @staticmethod
    def substitution(smile):

        # Returns where all the given benzene compound is substituted. If the given smile is not aromatic it simply returns.
        # The function returns a dict of the positions of various functional groups.

        # Example: c1cccc(Cl)c1Cl - Ortho dichlorobenzene.

        # Issues: So this code straight up wont work when the compound is some weird shiz with random stuff outside benzene or whatnot. This only works if
        # it is a simple benzene ring with the stuff inside validGroups as electrophiles.

        # Also, this does not work when the functional/electrophilic group is of the form
        # Brc1ccccc1
        # where the functional group Br appears before the benzene without brackets. I can fix this later as this is an edge case. For now, call it as
        # c1(Br)ccccc1

        if not Reactions.isAromatic(smile):
            # If the compound is not aromatic, we do not talk about its substitution.
            print("Given compund is not aromatic")
            return

        electrophilePosn = {}  # This dict will contain the posn of every electophile in the ring. It maps the carbon ring to position.
        # I start counting the carbon ring from the leftmost carbon atom to rightmost.

        # substring = "" # This is a string that keeps track of the current sub part of the given compund like a carbon atom/bromine etc.
        if not smile.startswith('c1'):
            startOfRing = smile.index('c1')
            electrophilePosn[1] = smile[0:startOfRing]

        carbonIndex = 0
        substring = ""
        for end in range(smile.index('c1'), len(smile)):  # Iterating through the smile.
            substring += smile[end]  # Adding the letter to substring
            # print(substring)
            if substring == 'c':  # If it is a small carbon atom (Benzene), we just increase carbon index.
                carbonIndex += 1
                substring = ""
            elif substring == '1':  # This is just a special case when u have a c1. I just split it into a special case by straight up ignoring the 1.
                substring = ""

            elif substring[0] == '(' and substring[-1] == ')' and substring.count(
                    '(') == substring.count(
                ')'):  # Now, if substring is a particular electrophile, then at 'carbonIndex', we have found an electrophile.
                if carbonIndex == 0:  # Special case when the electrophile happens in the first carbon in benzene. Carbon index is 0 in dat case so gotta change it.
                    electrophilePosn[carbonIndex + 1] = substring[1:-1]
                else:
                    electrophilePosn[carbonIndex] = substring[1:-1]

                substring = ""

        if substring != "":
            electrophilePosn[6] = substring

        return electrophilePosn

    @staticmethod
    def patternMatchElectrophile(smile):
        """
        Given a functional group, this function tries to pattern match it to a given direction for substitutiton.
        It returns a tuple where first index is the general form of the functional group and second index is the direction.
        Example, the group CCCCCC returns (R, O/P) as it is an alkyl group that directs to O/P.

        """
        # Alkyl groups
        if re.search("^C+$", smile):
            return ("R", "O/P")

        # Acetate. Not sure how to pattern match cuz idek how it looks like in smiles properly
        # if re.search("", smile):
        #   return ("", -1)

        # Ketones. Again not sure
        if re.search("^C*\(O\)C*$", smile):
            return ("COR", "M")

        return ("", -1)  # Default nothing matched

    @staticmethod
    def elecSub(smile, group):
        """
        Given a smile, this function will substitute it with the given group.
        """

        direction = {"C": "O/P", "O": "O/P", "N": "O/P", "F": "O/P", "Cl": "O/P", "Br/I": "O/P",
                     "[N+](=O)[O-]": "M", "C#N": "M", "C=O": "M", "O=C(O)": "M"}
        reactivity = []

        attachedGroups = Reactions.substitution(smile)

        if attachedGroups == {}:
            print("As only a single group is being added, we do not care about direction")
            return [group + smile]

        if len(attachedGroups) == 1:
            # Only a single func group attached
            item = attachedGroups.popitem()
            grp = item[1]
            posn = item[0]
            directionToSub = direction.get(grp, -1)  # Try to get the direction of the group

            groupAndDirectionFromFunc = None
            if directionToSub == -1:  # This means that the group is not there in the direction dictionary. So I will try some pattern matching to see
                # if there is any grp that fits.
                groupAndDirectionFromFunc = Reactions.patternMatchElectrophile(grp)
                directionToSub = groupAndDirectionFromFunc[1]

            # if groupAndDirectionFromFunc[0] == "": # Nothing is matched, then idk what to do lol.
            #   return -1 # Error

            if directionToSub == "O/P":
                positionToSub1 = ((posn) + 1) % 6  # Ortho position
                positionToSub2 = ((posn) + 3) % 6  # Para position

                noOfbenzylCarbon = 0  # Counts the number of benzyl carbon
                res = []  # Result storing the compounds produced

                for i in range(len(
                        smile)):  # Here, you are going through the smile and counting the no of carbons in the ring. When u find the position
                    # where u need to insert (After a given no of carbon atoms in benzene ring), you insert it.
                    if smile[i] == 'c':
                        noOfbenzylCarbon += 1

                    if noOfbenzylCarbon == positionToSub1:
                        index = i
                        if smile[i + 1] == '1':
                            index += 1

                        temp = "(" + group + ")"  # Adding brackets around group

                        res.append(smile[0: index + 1] + temp + smile[index + 1:])
                        break

                noOfbenzylCarbon = 0

                for i in range(len(smile)):
                    if smile[i] == 'c':
                        noOfbenzylCarbon += 1

                    if noOfbenzylCarbon == positionToSub2:
                        index = i
                        if smile[i + 1] == '1':
                            index += 1
                        temp = "(" + group + ")"

                        res.append(smile[0: index + 1] + temp + smile[index + 1:])
                        break

                return res

            elif directionToSub == "M":
                positionToSub1 = ((posn) + 2) % 6  # Meta position

                noOfbenzylCarbon = 0
                res = []

                for i in range(len(smile)):
                    if smile[i] == 'c':
                        noOfbenzylCarbon += 1

                    if noOfbenzylCarbon == positionToSub1:
                        index = i
                        if smile[i + 1] == '1':
                            index += 1

                        temp = "(" + group + ")"

                        res.append(smile[0: index + 1] + temp + smile[index + 1:])

                return res

        if len(attachedGroups) == 2:  # Two attached groups
            item1 = attachedGroups.popitem()
            item2 = attachedGroups.popitem()

    @staticmethod
    def chlorination(smile):
        return Reactions.elecSub(smile, "Cl")

    @staticmethod
    def florination(smile):
        return Reactions.elecSub(smile, "F")


""" 
This class will attempt to classify all possible reactions that can take place. The aim is to make the reactions 
as general as possible such that they are non-biased toward recognizing functional groups only in certain obvious 
circumstances and performing reactions upon them. I expect this will be a long-ass block of code.

This will contain the following reactions concerning alkenes (refer to the textbook for any clarification of concepts):

-Markovnikov addition of hydrogen halides to alkenes
-additions of halogens, X2, to alkenes
-halohydrin formation across the alkene double-bond: X2 + H2O add with Markovnikov-regiochemistry and anti-stereochemistry
-addition of water by oxymercuration-demercuration (Markovnikov regiochemistry)
-addition of water by hydroboration-oxidation (non-Markovnikov regiochemistry)
-catalytic hydrogenation
-epoxidation with a peroxyacid
-hydroxylation with osmium tetroxide
-dichlorocarbene addition to yield cyclopropanes
-Simmons-Smith reaction to yield cyclopropanes
-hydroxylation by acid-catalyzed epoxide hydrolysis
-oxidative cleavage of alkenes by reaction with ozone followed by zinc in acetic acid
-oxidative cleavage of alkenes by reaction with potassium permanganate in acidic solution
-cleavage of 1,2-diols

"""


class AlkeneReactions:
    """
    The functions HBr_addition, HCl_addition, and HI_addition are identical, except for the halogens of course.
    """

    @staticmethod
    def HBr_addition(smile):
        smile_products = [
            smile]  # this list will contain all possible products, accounting for stereochemistry
        #draw(smile)

        smile1 = smile_products[0]

        while FunctionalGroups.alkene(smile1) == True:

            # We must consider regioselectivity for hydrogen halide addition to an alkene. The halide will add to the more
            # substituted alkene carbon and the hydrogen will add to the less-substituted one. In the case where both alkene
            # carbons have equal number of carbon substituents, two possible products can form since the halide can bond
            # to either carbon.

            # one alkene carbon has two carbons singly-bonded to it; the other alkene carbon is not bonded to any carbons
            rxn1 = AllChem.ReactionFromSmarts(
                '[C:0][C:1]([C:4])=[CH2:2].[Br:3]>>[C:0][C:1]([C:4])(-[Br:3])[CH2:2]',
                useSmiles=True)

            # This rxn format takes in a SMARTS code that can find, generally, alkenes of a certain degree of substitution.
            # ex. In rxn1, [CH0:1]([C:4])=[CH2:2] indicates that a disubstituted alkene carbon, "C1" with 0 H bonded to it,
            # is doubly-bonded to a nonsubstituted alkene carbon, "C2" with 2 H's bonded to it. The "." separates two
            # unbonded compounds, in this case, the alkene and the bromine. The reaction symbol ">>" then indicates reaction,
            # from which results the halogen addition product.

            ps1 = rxn1.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H][Br]")))
            # This ps format gives the SMILES output of the reaction product between the SMILES input and HBr.

            # one alkene carbon is singly-bonded to two carbons, the other is singly-bonded to one carbon
            rxn2 = AllChem.ReactionFromSmarts(
                '[C:6][C:1]([C:4])=[CH:2]([C:5]).[Br:3]>>[C:6][C:1]([C:4])(-[Br:3])[CH:2]([C:5])',
                useSmiles=True)
            ps2 = rxn2.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H][Br]")))

            # both alkene carbons are each singly-bonded to two carbons
            rxn3 = AllChem.ReactionFromSmarts(
                '[C:0][C:1]([C:2])=[C:3]([C:4])[C:5].[Br:6]>>[C:0][C:1]([C:2])(-[Br:6])[C:3]([C:4])[C:5]',
                useSmiles=True)
            ps3 = rxn3.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H][Br]")))

            # one alkene carbon is singly-bonded to one carbon, the other is not bonded to any carbons
            rxn4 = AllChem.ReactionFromSmarts(
                '[C:4][CH:1]=[CH2:2].[Br:3]>>[C:4][CH:1](-[Br:3])-[CH2:2]', useSmiles=True)
            ps4 = rxn4.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H][Br]")))

            # both alkene carbons each singly-bonded to one carbon
            rxn5 = AllChem.ReactionFromSmarts(
                '[C:4][CH:1]=[CH:2][C:5].[Br:3]>>[C:4][CH:1](-[Br:3])[CH:2][C:5]', useSmiles=True)
            ps5 = rxn5.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H][Br]")))

            # both alkene carbons are not bonded to any carbons
            rxn6 = AllChem.ReactionFromSmarts('[CH2:1]=[CH2:2].[Br:3]>>[CH2:1](-[Br:3])[CH3:2]',
                                              useSmiles=True)
            ps6 = rxn6.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H][Br]")))

            m1 = Chem.MolFromSmiles(smile1)
            patt1 = Chem.MolFromSmarts("[C][C]([C])=[CH2]")
            r1 = m1.HasSubstructMatch(patt1)

            m2 = Chem.MolFromSmiles(smile1)
            patt2 = Chem.MolFromSmarts("[C][C]([C])=[CH][C]")
            r2 = m2.HasSubstructMatch(patt2)

            m3 = Chem.MolFromSmiles(smile1)
            patt3 = Chem.MolFromSmarts("[C][C]([C])=[C]([C])[C]")
            r3 = m3.HasSubstructMatch(patt3)

            m4 = Chem.MolFromSmiles(smile1)
            patt4 = Chem.MolFromSmarts("[C][CH]=[CH2]")
            r4 = m4.HasSubstructMatch(patt4)

            m5 = Chem.MolFromSmiles(smile1)
            patt5 = Chem.MolFromSmarts("[C][CH]=[CH][C]")
            r5 = m5.HasSubstructMatch(patt5)

            m6 = Chem.MolFromSmiles(smile1)
            patt6 = Chem.MolFromSmarts("[CH2]=[CH2]")
            r6 = m6.HasSubstructMatch(patt6)

            alkene_H_halide_reactions_list = [rxn1, rxn2, rxn3, rxn4, rxn5, rxn6]
            H_halideaddition_products_list = [ps1, ps2, ps3, ps4, ps5, ps6]
            alkene_regiochemistry_list = [r1, r2, r3, r4, r5, r6]
            num_possible_products_list = [1, 1, 2, 1, 2, 1]

            i = 0
            for i in range(6):
                if alkene_regiochemistry_list[i] == True:
                    j = 0
                    for j in range(num_possible_products_list[i]):
                        smile_products.append(
                            Chem.MolToSmiles(H_halideaddition_products_list[i][j][0]))
                        j += 1
                    del (smile_products[0])
                    smile1 = smile_products[0]

                    i = 0
                    break

                i += 1

        k = 0
        # for k in range(len(smile_products)):
        #     draw(smile_products[k])
        #     k += 1

        return smile_products

    @staticmethod
    def HCl_addition(smile):
        smile_products = [smile]
        #draw(smile)

        smile1 = smile_products[0]

        while FunctionalGroups.alkene(smile1) == True:

            # one alkene carbon has two carbons singly-bonded to it; the other alkene carbon is not bonded to any carbons
            rxn1 = AllChem.ReactionFromSmarts(
                '[C:0][C:1]([C:4])=[CH2:2].[Cl:3]>>[C:0][C:1]([C:4])(-[Cl:3])[CH2:2]',
                useSmiles=True)
            ps1 = rxn1.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H][Cl]")))

            # one alkene carbon is singly-bonded to two carbons, the other is singly-bonded to one carbon
            rxn2 = AllChem.ReactionFromSmarts(
                '[C:6][C:1]([C:4])=[CH:2]([C:5]).[Cl:3]>>[C:6][C:1]([C:4])(-[Cl:3])[CH:2]([C:5])',
                useSmiles=True)
            ps2 = rxn2.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H][Cl]")))

            # both alkene carbons are each singly-bonded to two carbons
            rxn3 = AllChem.ReactionFromSmarts(
                '[C:0][C:1]([C:2])=[C:3]([C:4])[C:5].[Cl:6]>>[C:0][C:1]([C:2])(-[Cl:6])[C:3]([C:4])[C:5]',
                useSmiles=True)
            ps3 = rxn3.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H][Cl]")))

            # one alkene carbon is singly-bonded to one carbon, the other is not bonded to any carbons
            rxn4 = AllChem.ReactionFromSmarts(
                '[C:4][CH:1]=[CH2:2].[Cl:3]>>[C:4][CH:1](-[Cl:3])-[CH2:2]', useSmiles=True)
            ps4 = rxn4.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H][Cl]")))

            # both alkene carbons each singly-bonded to one carbon
            rxn5 = AllChem.ReactionFromSmarts(
                '[C:4][CH:1]=[CH:2][C:5].[Cl:3]>>[C:4][CH:1](-[Cl:3])[CH:2][C:5]', useSmiles=True)
            ps5 = rxn5.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H][Cl]")))

            # both alkene carbons are not bonded to any carbons
            rxn6 = AllChem.ReactionFromSmarts('[CH2:1]=[CH2:2].[Cl:3]>>[CH2:1](-[Cl:3])[CH3:2]',
                                              useSmiles=True)
            ps6 = rxn6.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H][Cl]")))

            m1 = Chem.MolFromSmiles(smile1)
            patt1 = Chem.MolFromSmarts("[C][C]([C])=[CH2]")
            r1 = m1.HasSubstructMatch(patt1)

            m2 = Chem.MolFromSmiles(smile1)
            patt2 = Chem.MolFromSmarts("[C][C]([C])=[CH][C]")
            r2 = m2.HasSubstructMatch(patt2)

            m3 = Chem.MolFromSmiles(smile1)
            patt3 = Chem.MolFromSmarts("[C][C]([C])=[C]([C])[C]")
            r3 = m3.HasSubstructMatch(patt3)

            m4 = Chem.MolFromSmiles(smile1)
            patt4 = Chem.MolFromSmarts("[C][CH]=[CH2]")
            r4 = m4.HasSubstructMatch(patt4)

            m5 = Chem.MolFromSmiles(smile1)
            patt5 = Chem.MolFromSmarts("[C][CH]=[CH][C]")
            r5 = m5.HasSubstructMatch(patt5)

            m6 = Chem.MolFromSmiles(smile1)
            patt6 = Chem.MolFromSmarts("[CH2]=[CH2]")
            r6 = m6.HasSubstructMatch(patt6)

            alkene_H_halide_reactions_list = [rxn1, rxn2, rxn3, rxn4, rxn5, rxn6]
            H_halideaddition_products_list = [ps1, ps2, ps3, ps4, ps5, ps6]
            alkene_regiochemistry_list = [r1, r2, r3, r4, r5, r6]
            num_possible_products_list = [1, 1, 2, 1, 2, 1]

            i = 0
            for i in range(6):
                if alkene_regiochemistry_list[i] == True:
                    j = 0
                    for j in range(num_possible_products_list[i]):
                        smile_products.append(
                            Chem.MolToSmiles(H_halideaddition_products_list[i][j][0]))
                        j += 1
                    del (smile_products[0])
                    smile1 = smile_products[0]

                    i = 0
                    break

                i += 1

        k = 0
        # for k in range(len(smile_products)):
        #     draw(smile_products[k])
        #     k += 1

        return smile_products

    @staticmethod
    def HI_addition(smile):
        smile_products = [smile]
        #draw(smile)

        smile1 = smile_products[0]

        while FunctionalGroups.alkene(smile1) == True:

            # one alkene carbon has two carbons singly-bonded to it; the other alkene carbon is not bonded to any carbons
            rxn1 = AllChem.ReactionFromSmarts(
                '[C:0][C:1]([C:4])=[CH2:2].[I:3]>>[C:0][C:1]([C:4])(-[I:3])[CH2:2]', useSmiles=True)
            ps1 = rxn1.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H][I]")))

            # one alkene carbon is singly-bonded to two carbons, the other is singly-bonded to one carbon
            rxn2 = AllChem.ReactionFromSmarts(
                '[C:6][C:1]([C:4])=[CH:2]([C:5]).[I:3]>>[C:6][C:1]([C:4])(-[I:3])[CH:2]([C:5])',
                useSmiles=True)
            ps2 = rxn2.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H][I]")))

            # both alkene carbons are each singly-bonded to two carbons
            rxn3 = AllChem.ReactionFromSmarts(
                '[C:0][C:1]([C:2])=[C:3]([C:4])[C:5].[I:6]>>[C:0][C:1]([C:2])(-[I:6])[C:3]([C:4])[C:5]',
                useSmiles=True)
            ps3 = rxn3.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H][I]")))

            # one alkene carbon is singly-bonded to one carbon, the other is not bonded to any carbons
            rxn4 = AllChem.ReactionFromSmarts(
                '[C:4][CH:1]=[CH2:2].[I:3]>>[C:4][CH:1](-[I:3])-[CH2:2]', useSmiles=True)
            ps4 = rxn4.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H][I]")))

            # both alkene carbons each singly-bonded to one carbon
            rxn5 = AllChem.ReactionFromSmarts(
                '[C:4][CH:1]=[CH:2][C:5].[I:3]>>[C:4][CH:1](-[I:3])[CH:2][C:5]', useSmiles=True)
            ps5 = rxn5.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H][I]")))

            # both alkene carbons are not bonded to any carbons
            rxn6 = AllChem.ReactionFromSmarts('[CH2:1]=[CH2:2].[I:3]>>[CH2:1](-[I:3])[CH3:2]',
                                              useSmiles=True)
            ps6 = rxn6.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H][I]")))

            m1 = Chem.MolFromSmiles(smile1)
            patt1 = Chem.MolFromSmarts("[C][C]([C])=[CH2]")
            r1 = m1.HasSubstructMatch(patt1)

            m2 = Chem.MolFromSmiles(smile1)
            patt2 = Chem.MolFromSmarts("[C][C]([C])=[CH][C]")
            r2 = m2.HasSubstructMatch(patt2)

            m3 = Chem.MolFromSmiles(smile1)
            patt3 = Chem.MolFromSmarts("[C][C]([C])=[C]([C])[C]")
            r3 = m3.HasSubstructMatch(patt3)

            m4 = Chem.MolFromSmiles(smile1)
            patt4 = Chem.MolFromSmarts("[C][CH]=[CH2]")
            r4 = m4.HasSubstructMatch(patt4)

            m5 = Chem.MolFromSmiles(smile1)
            patt5 = Chem.MolFromSmarts("[C][CH]=[CH][C]")
            r5 = m5.HasSubstructMatch(patt5)

            m6 = Chem.MolFromSmiles(smile1)
            patt6 = Chem.MolFromSmarts("[CH2]=[CH2]")
            r6 = m6.HasSubstructMatch(patt6)

            alkene_H_halide_reactions_list = [rxn1, rxn2, rxn3, rxn4, rxn5, rxn6]
            H_halideaddition_products_list = [ps1, ps2, ps3, ps4, ps5, ps6]
            alkene_regiochemistry_list = [r1, r2, r3, r4, r5, r6]
            num_possible_products_list = [1, 1, 2, 1, 2, 1]

            i = 0
            for i in range(6):
                if alkene_regiochemistry_list[i] == True:
                    j = 0
                    for j in range(num_possible_products_list[i]):
                        smile_products.append(
                            Chem.MolToSmiles(H_halideaddition_products_list[i][j][0]))
                        j += 1
                    del (smile_products[0])
                    smile1 = smile_products[0]

                    i = 0
                    break

                i += 1

        k = 0
        # for k in range(len(smile_products)):
        #     draw(smile_products[k])
        #     k += 1

        return smile_products

    @staticmethod
    def Br2_addition(smile):
        smile_products = [smile]
        #draw(smile)

        smile1 = smile_products[0]

        while FunctionalGroups.alkene(smile1) == True:
            rxn1 = AllChem.ReactionFromSmarts(
                '[C:1]=[C:2].[Br:3][Br:4]>>[C:1](-[Br:3])[C:2](-[Br:4])', useSmiles=True)
            ps1 = rxn1.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("BrBr")))
            smile_products.append(Chem.MolToSmiles(ps1[0][0]))
            del (smile_products[0])
            smile1 = smile_products[0]

        k = 0
        # for k in range(len(smile_products)):
        #     draw(smile_products[k])
        #     k += 1

        return smile_products

    @staticmethod
    def Cl2_addition(smile):
        smile_products = [smile]
        # draw(smile)

        smile1 = smile_products[0]

        while FunctionalGroups.alkene(smile1) == True:
            rxn1 = AllChem.ReactionFromSmarts(
                '[C:1]=[C:2].[Cl:3][Cl:4]>>[C:1](-[Cl:3])[C:2](-[Cl:4])', useSmiles=True)
            ps1 = rxn1.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("ClCl")))
            smile_products.append(Chem.MolToSmiles(ps1[0][0]))
            del (smile_products[0])
            smile1 = smile_products[0]

        k = 0
        for k in range(len(smile_products)):
            # draw(smile_products[k])
            k += 1
        return smile_products

    @staticmethod
    def bromohydrin_formation(smile):  # the "-OH" will bond to the more-substituted carbon
        smile_products = [smile]
        # draw(smile)

        smile1 = smile_products[0]

        while FunctionalGroups.alkene(smile1) == True:

            # one alkene carbon has two carbons singly-bonded to it; the other alkene carbon is not bonded to any carbons
            rxn1 = AllChem.ReactionFromSmarts(
                '[C:0][C:1]([C:4])=[CH2:2].[O:8][Br:9]>>[C:0][C:1]([C:4])(-[O:8])[CH2:2](-[Br:9])',
                useSmiles=True)
            ps1 = rxn1.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("BrO")))

            # one alkene carbon is singly-bonded to two carbons, the other is singly-bonded to one carbon
            rxn2 = AllChem.ReactionFromSmarts(
                '[C:6][C:1]([C:4])=[CH:2]([C:5]).[O:8][Br:9]>>[C:6][C:1]([C:4])(-[O:8])[CH:2](-[Br:9])([C:5])',
                useSmiles=True)
            ps2 = rxn2.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("BrO")))

            # both alkene carbons are each singly-bonded to two carbons
            rxn3 = AllChem.ReactionFromSmarts(
                '[C:0][C:1]([C:2])=[C:3]([C:4])[C:5].[O:8][Br:9]>>[C:0][C:1]([C:2])(-[O:8])[C:3](-[Br:9])([C:4])[C:5]',
                useSmiles=True)
            ps3 = rxn3.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("BrO")))

            # one alkene carbon is singly-bonded to one carbon, the other is not bonded to any carbons
            rxn4 = AllChem.ReactionFromSmarts(
                '[C:4][CH:1]=[CH2:2].[O:8][Br:9]>>[C:4][CH:1](-[O:8])-[CH2:2](-[Br:9])',
                useSmiles=True)
            ps4 = rxn4.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("BrO")))

            # both alkene carbons each singly-bonded to one carbon
            rxn5 = AllChem.ReactionFromSmarts(
                '[C:4][CH:1]=[CH:2][C:5].[O:8][Br:9]>>[C:4][CH:1](-[O:8])[CH:2](-[Br:9])[C:5]',
                useSmiles=True)
            ps5 = rxn5.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("BrO")))

            # both alkene carbons are not bonded to any carbons
            rxn6 = AllChem.ReactionFromSmarts(
                '[CH2:1]=[CH2:2].[O:8][Br:9]>>[CH2:1](-[O:8])[CH3:2](-[Br:9])', useSmiles=True)
            ps6 = rxn6.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("BrO")))

            m1 = Chem.MolFromSmiles(smile1)
            patt1 = Chem.MolFromSmarts("[C][C]([C])=[CH2]")
            r1 = m1.HasSubstructMatch(patt1)

            m2 = Chem.MolFromSmiles(smile1)
            patt2 = Chem.MolFromSmarts("[C][C]([C])=[CH][C]")
            r2 = m2.HasSubstructMatch(patt2)

            m3 = Chem.MolFromSmiles(smile1)
            patt3 = Chem.MolFromSmarts("[C][C]([C])=[C]([C])[C]")
            r3 = m3.HasSubstructMatch(patt3)

            m4 = Chem.MolFromSmiles(smile1)
            patt4 = Chem.MolFromSmarts("[C][CH]=[CH2]")
            r4 = m4.HasSubstructMatch(patt4)

            m5 = Chem.MolFromSmiles(smile1)
            patt5 = Chem.MolFromSmarts("[C][CH]=[CH][C]")
            r5 = m5.HasSubstructMatch(patt5)

            m6 = Chem.MolFromSmiles(smile1)
            patt6 = Chem.MolFromSmarts("[CH2]=[CH2]")
            r6 = m6.HasSubstructMatch(patt6)

            alkene_H_halide_reactions_list = [rxn1, rxn2, rxn3, rxn4, rxn5, rxn6]
            H_halideaddition_products_list = [ps1, ps2, ps3, ps4, ps5, ps6]
            alkene_regiochemistry_list = [r1, r2, r3, r4, r5, r6]
            num_possible_products_list = [1, 1, 2, 1, 2, 1]

            i = 0
            for i in range(6):
                if alkene_regiochemistry_list[i] == True:
                    j = 0
                    for j in range(num_possible_products_list[i]):
                        smile_products.append(
                            Chem.MolToSmiles(H_halideaddition_products_list[i][j][0]))
                        j += 1
                    del (smile_products[0])
                    smile1 = smile_products[0]

                    i = 0
                    break

                i += 1

        k = 0
        # for k in range(len(smile_products)):
        #     draw(smile_products[k])
        #     k += 1

        return smile_products

    @staticmethod
    def chlorohydrin_formation(smile):
        smile_products = [smile]
        # draw(smile)

        smile1 = smile_products[0]

        while FunctionalGroups.alkene(smile1) == True:

            # one alkene carbon has two carbons singly-bonded to it; the other alkene carbon is not bonded to any carbons
            rxn1 = AllChem.ReactionFromSmarts(
                '[C:0][C:1]([C:4])=[CH2:2].[O:8][Cl:9]>>[C:0][C:1]([C:4])(-[O:8])[CH2:2](-[Cl:9])',
                useSmiles=True)
            ps1 = rxn1.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("ClO")))

            # one alkene carbon is singly-bonded to two carbons, the other is singly-bonded to one carbon
            rxn2 = AllChem.ReactionFromSmarts(
                '[C:6][C:1]([C:4])=[CH:2]([C:5]).[O:8][Cl:9]>>[C:6][C:1]([C:4])(-[O:8])[CH:2](-[Cl:9])([C:5])',
                useSmiles=True)
            ps2 = rxn2.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("ClO")))

            # both alkene carbons are each singly-bonded to two carbons
            rxn3 = AllChem.ReactionFromSmarts(
                '[C:0][C:1]([C:2])=[C:3]([C:4])[C:5].[O:8][Cl:9]>>[C:0][C:1]([C:2])(-[O:8])[C:3](-[Cl:9])([C:4])[C:5]',
                useSmiles=True)
            ps3 = rxn3.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("ClO")))

            # one alkene carbon is singly-bonded to one carbon, the other is not bonded to any carbons
            rxn4 = AllChem.ReactionFromSmarts(
                '[C:4][CH:1]=[CH2:2].[O:8][Cl:9]>>[C:4][CH:1](-[O:8])-[CH2:2](-[Cl:9])',
                useSmiles=True)
            ps4 = rxn4.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("ClO")))

            # both alkene carbons each singly-bonded to one carbon
            rxn5 = AllChem.ReactionFromSmarts(
                '[C:4][CH:1]=[CH:2][C:5].[O:8][Cl:9]>>[C:4][CH:1](-[O:8])[CH:2](-[Cl:9])[C:5]',
                useSmiles=True)
            ps5 = rxn5.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("ClO")))

            # both alkene carbons are not bonded to any carbons
            rxn6 = AllChem.ReactionFromSmarts(
                '[CH2:1]=[CH2:2].[O:8][Cl:9]>>[CH2:1](-[O:8])[CH3:2](-[Cl:9])', useSmiles=True)
            ps6 = rxn6.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("ClO")))

            m1 = Chem.MolFromSmiles(smile1)
            patt1 = Chem.MolFromSmarts("[C][C]([C])=[CH2]")
            r1 = m1.HasSubstructMatch(patt1)

            m2 = Chem.MolFromSmiles(smile1)
            patt2 = Chem.MolFromSmarts("[C][C]([C])=[CH][C]")
            r2 = m2.HasSubstructMatch(patt2)

            m3 = Chem.MolFromSmiles(smile1)
            patt3 = Chem.MolFromSmarts("[C][C]([C])=[C]([C])[C]")
            r3 = m3.HasSubstructMatch(patt3)

            m4 = Chem.MolFromSmiles(smile1)
            patt4 = Chem.MolFromSmarts("[C][CH]=[CH2]")
            r4 = m4.HasSubstructMatch(patt4)

            m5 = Chem.MolFromSmiles(smile1)
            patt5 = Chem.MolFromSmarts("[C][CH]=[CH][C]")
            r5 = m5.HasSubstructMatch(patt5)

            m6 = Chem.MolFromSmiles(smile1)
            patt6 = Chem.MolFromSmarts("[CH2]=[CH2]")
            r6 = m6.HasSubstructMatch(patt6)

            alkene_H_halide_reactions_list = [rxn1, rxn2, rxn3, rxn4, rxn5, rxn6]
            H_halideaddition_products_list = [ps1, ps2, ps3, ps4, ps5, ps6]
            alkene_regiochemistry_list = [r1, r2, r3, r4, r5, r6]
            num_possible_products_list = [1, 1, 2, 1, 2, 1]

            i = 0
            for i in range(6):
                if alkene_regiochemistry_list[i] == True:
                    j = 0
                    for j in range(num_possible_products_list[i]):
                        smile_products.append(
                            Chem.MolToSmiles(H_halideaddition_products_list[i][j][0]))
                        j += 1
                    del (smile_products[0])
                    smile1 = smile_products[0]

                    i = 0
                    break

                i += 1

        k = 0
        # for k in range(len(smile_products)):
        #     draw(smile_products[k])
        #     k += 1

        return smile_products

    @staticmethod
    def iodohydrin_formation(smile):
        smile_products = [smile]
        # draw(smile)

        smile1 = smile_products[0]

        while FunctionalGroups.alkene(smile1) == True:

            # one alkene carbon has two carbons singly-bonded to it; the other alkene carbon is not bonded to any carbons
            rxn1 = AllChem.ReactionFromSmarts(
                '[C:0][C:1]([C:4])=[CH2:2].[O:8][I:9]>>[C:0][C:1]([C:4])(-[O:8])[CH2:2](-[I:9])',
                useSmiles=True)
            ps1 = rxn1.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("IO")))

            # one alkene carbon is singly-bonded to two carbons, the other is singly-bonded to one carbon
            rxn2 = AllChem.ReactionFromSmarts(
                '[C:6][C:1]([C:4])=[CH:2]([C:5]).[O:8][I:9]>>[C:6][C:1]([C:4])(-[O:8])[CH:2](-[I:9])([C:5])',
                useSmiles=True)
            ps2 = rxn2.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("IO")))

            # both alkene carbons are each singly-bonded to two carbons
            rxn3 = AllChem.ReactionFromSmarts(
                '[C:0][C:1]([C:2])=[C:3]([C:4])[C:5].[O:8][I:9]>>[C:0][C:1]([C:2])(-[O:8])[C:3](-[I:9])([C:4])[C:5]',
                useSmiles=True)
            ps3 = rxn3.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("IO")))

            # one alkene carbon is singly-bonded to one carbon, the other is not bonded to any carbons
            rxn4 = AllChem.ReactionFromSmarts(
                '[C:4][CH:1]=[CH2:2].[O:8][I:9]>>[C:4][CH:1](-[O:8])-[CH2:2](-[I:9])',
                useSmiles=True)
            ps4 = rxn4.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("IO")))

            # both alkene carbons each singly-bonded to one carbon
            rxn5 = AllChem.ReactionFromSmarts(
                '[C:4][CH:1]=[CH:2][C:5].[O:8][I:9]>>[C:4][CH:1](-[O:8])[CH:2](-[I:9])[C:5]',
                useSmiles=True)
            ps5 = rxn5.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("IO")))

            # both alkene carbons are not bonded to any carbons
            rxn6 = AllChem.ReactionFromSmarts(
                '[CH2:1]=[CH2:2].[O:8][I:9]>>[CH2:1](-[O:8])[CH3:2](-[I:9])', useSmiles=True)
            ps6 = rxn6.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("IO")))

            m1 = Chem.MolFromSmiles(smile1)
            patt1 = Chem.MolFromSmarts("[C][C]([C])=[CH2]")
            r1 = m1.HasSubstructMatch(patt1)

            m2 = Chem.MolFromSmiles(smile1)
            patt2 = Chem.MolFromSmarts("[C][C]([C])=[CH][C]")
            r2 = m2.HasSubstructMatch(patt2)

            m3 = Chem.MolFromSmiles(smile1)
            patt3 = Chem.MolFromSmarts("[C][C]([C])=[C]([C])[C]")
            r3 = m3.HasSubstructMatch(patt3)

            m4 = Chem.MolFromSmiles(smile1)
            patt4 = Chem.MolFromSmarts("[C][CH]=[CH2]")
            r4 = m4.HasSubstructMatch(patt4)

            m5 = Chem.MolFromSmiles(smile1)
            patt5 = Chem.MolFromSmarts("[C][CH]=[CH][C]")
            r5 = m5.HasSubstructMatch(patt5)

            m6 = Chem.MolFromSmiles(smile1)
            patt6 = Chem.MolFromSmarts("[CH2]=[CH2]")
            r6 = m6.HasSubstructMatch(patt6)

            alkene_H_halide_reactions_list = [rxn1, rxn2, rxn3, rxn4, rxn5, rxn6]
            H_halideaddition_products_list = [ps1, ps2, ps3, ps4, ps5, ps6]
            alkene_regiochemistry_list = [r1, r2, r3, r4, r5, r6]
            num_possible_products_list = [1, 1, 2, 1, 2, 1]

            i = 0
            for i in range(6):
                if alkene_regiochemistry_list[i] == True:
                    j = 0
                    for j in range(num_possible_products_list[i]):
                        smile_products.append(
                            Chem.MolToSmiles(H_halideaddition_products_list[i][j][0]))
                        j += 1
                    del (smile_products[0])
                    smile1 = smile_products[0]

                    i = 0
                    break

                i += 1

        k = 0
        # for k in range(len(smile_products)):
        #     draw(smile_products[k])
        #     k += 1

        return smile_products

    @staticmethod
    def Markovnikov_water_addition(smile):
        smile_products = [
            smile]  # following Markovnikov's rule, the "-OH" will ad to the more-substituted carbon
        # draw(smile)

        smile1 = smile_products[0]

        while FunctionalGroups.alkene(smile1) == True:
            # one alkene carbon has two carbons singly-bonded to it; the other alkene carbon is not bonded to any carbons
            rxn1 = AllChem.ReactionFromSmarts(
                '[C:0][C:1]([C:4])=[CH2:2].[O:3]>>[C:0][C:1]([C:4])(-[O:3])[CH2:2]', useSmiles=True)
            ps1 = rxn1.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H]O[H]")))

            # one alkene carbon is singly-bonded to two carbons, the other is singly-bonded to one carbon
            rxn2 = AllChem.ReactionFromSmarts(
                '[C:6][C:1]([C:4])=[CH:2]([C:5]).[O:3]>>[C:6][C:1]([C:4])(-[O:3])[CH:2]([C:5])',
                useSmiles=True)
            ps2 = rxn2.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H]O[H]")))

            # both alkene carbons are each singly-bonded to two carbons
            rxn3 = AllChem.ReactionFromSmarts(
                '[C:0][C:1]([C:2])=[C:3]([C:4])[C:5].[O:6]>>[C:0][C:1]([C:2])(-[O:6])[C:3]([C:4])[C:5]',
                useSmiles=True)
            ps3 = rxn3.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H]O[H]")))

            # one alkene carbon is singly-bonded to one carbon, the other is not bonded to any carbons
            rxn4 = AllChem.ReactionFromSmarts(
                '[C:4][CH:1]=[CH2:2].[O:3]>>[C:4][CH:1](-[O:3])-[CH2:2]', useSmiles=True)
            ps4 = rxn4.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H]O[H]")))

            # both alkene carbons each singly-bonded to one carbon
            rxn5 = AllChem.ReactionFromSmarts(
                '[C:4][CH:1]=[CH:2][C:5].[O:3]>>[C:4][CH:1](-[O:3])[CH:2][C:5]', useSmiles=True)
            ps5 = rxn5.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H]O[H]")))

            # both alkene carbons are not bonded to any carbons
            rxn6 = AllChem.ReactionFromSmarts('[CH2:1]=[CH2:2].[O:3]>>[CH2:1](-[O:3])[CH3:2]',
                                              useSmiles=True)
            ps6 = rxn6.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H]O[H]")))

            m1 = Chem.MolFromSmiles(smile1)
            patt1 = Chem.MolFromSmarts("[C][C]([C])=[CH2]")
            r1 = m1.HasSubstructMatch(patt1)

            m2 = Chem.MolFromSmiles(smile1)
            patt2 = Chem.MolFromSmarts("[C][C]([C])=[CH][C]")
            r2 = m2.HasSubstructMatch(patt2)

            m3 = Chem.MolFromSmiles(smile1)
            patt3 = Chem.MolFromSmarts("[C][C]([C])=[C]([C])[C]")
            r3 = m3.HasSubstructMatch(patt3)

            m4 = Chem.MolFromSmiles(smile1)
            patt4 = Chem.MolFromSmarts("[C][CH]=[CH2]")
            r4 = m4.HasSubstructMatch(patt4)

            m5 = Chem.MolFromSmiles(smile1)
            patt5 = Chem.MolFromSmarts("[C][CH]=[CH][C]")
            r5 = m5.HasSubstructMatch(patt5)

            m6 = Chem.MolFromSmiles(smile1)
            patt6 = Chem.MolFromSmarts("[CH2]=[CH2]")
            r6 = m6.HasSubstructMatch(patt6)

            alkene_H_halide_reactions_list = [rxn1, rxn2, rxn3, rxn4, rxn5, rxn6]
            H_halideaddition_products_list = [ps1, ps2, ps3, ps4, ps5, ps6]
            alkene_regiochemistry_list = [r1, r2, r3, r4, r5, r6]
            num_possible_products_list = [1, 1, 2, 1, 2, 1]

            i = 0
            for i in range(6):
                if alkene_regiochemistry_list[i] == True:
                    j = 0
                    for j in range(num_possible_products_list[i]):
                        smile_products.append(
                            Chem.MolToSmiles(H_halideaddition_products_list[i][j][0]))
                        j += 1
                    del (smile_products[0])
                    smile1 = smile_products[0]

                    i = 0
                    break

                i += 1

        k = 0
        # for k in range(len(smile_products)):
        #     draw(smile_products[k])
        #     k += 1

        return smile_products

    @staticmethod
    def nonMarkovnikov_water_addition(smile):
        smile_products = [smile]  # the "-OH" will add to the less-substituted carbon
        # draw(smile)

        smile1 = smile_products[0]

        while FunctionalGroups.alkene(smile1) == True:
            # one alkene carbon has two carbons singly-bonded to it; the other alkene carbon is not bonded to any carbons
            rxn1 = AllChem.ReactionFromSmarts(
                '[C:0][C:1]([C:4])=[CH2:2].[O:3]>>[C:0][C:1]([C:4])[CH2:2](-[O:3])', useSmiles=True)
            ps1 = rxn1.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H]O[H]")))

            # one alkene carbon is singly-bonded to two carbons, the other is singly-bonded to one carbon
            rxn2 = AllChem.ReactionFromSmarts(
                '[C:6][C:1]([C:4])=[CH:2]([C:5]).[O:3]>>[C:6][C:1]([C:4])[CH:2](-[O:3])([C:5])',
                useSmiles=True)
            ps2 = rxn2.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H]O[H]")))

            # both alkene carbons are each singly-bonded to two carbons
            rxn3 = AllChem.ReactionFromSmarts(
                '[C:0][C:1]([C:2])=[C:3]([C:4])[C:5].[O:6]>>[C:0][C:1]([C:2])[C:3](-[O:6])([C:4])[C:5]',
                useSmiles=True)
            ps3 = rxn3.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H]O[H]")))

            # one alkene carbon is singly-bonded to one carbon, the other is not bonded to any carbons
            rxn4 = AllChem.ReactionFromSmarts(
                '[C:4][CH:1]=[CH2:2].[O:3]>>[C:4][CH:1]-[CH2:2](-[O:3])', useSmiles=True)
            ps4 = rxn4.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H]O[H]")))

            # both alkene carbons each singly-bonded to one carbon
            rxn5 = AllChem.ReactionFromSmarts(
                '[C:4][CH:1]=[CH:2][C:5].[O:3]>>[C:4][CH:1][CH:2](-[O:3])[C:5]', useSmiles=True)
            ps5 = rxn5.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H]O[H]")))

            # both alkene carbons are not bonded to any carbons
            rxn6 = AllChem.ReactionFromSmarts('[CH2:1]=[CH2:2].[O:3]>>[CH2:1](-[O:3])[CH3:2]',
                                              useSmiles=True)
            ps6 = rxn6.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H]O[H]")))

            m1 = Chem.MolFromSmiles(smile1)
            patt1 = Chem.MolFromSmarts("[C][C]([C])=[CH2]")
            r1 = m1.HasSubstructMatch(patt1)

            m2 = Chem.MolFromSmiles(smile1)
            patt2 = Chem.MolFromSmarts("[C][C]([C])=[CH][C]")
            r2 = m2.HasSubstructMatch(patt2)

            m3 = Chem.MolFromSmiles(smile1)
            patt3 = Chem.MolFromSmarts("[C][C]([C])=[C]([C])[C]")
            r3 = m3.HasSubstructMatch(patt3)

            m4 = Chem.MolFromSmiles(smile1)
            patt4 = Chem.MolFromSmarts("[C][CH]=[CH2]")
            r4 = m4.HasSubstructMatch(patt4)

            m5 = Chem.MolFromSmiles(smile1)
            patt5 = Chem.MolFromSmarts("[C][CH]=[CH][C]")
            r5 = m5.HasSubstructMatch(patt5)

            m6 = Chem.MolFromSmiles(smile1)
            patt6 = Chem.MolFromSmarts("[CH2]=[CH2]")
            r6 = m6.HasSubstructMatch(patt6)

            alkene_H_halide_reactions_list = [rxn1, rxn2, rxn3, rxn4, rxn5, rxn6]
            H_halideaddition_products_list = [ps1, ps2, ps3, ps4, ps5, ps6]
            alkene_regiochemistry_list = [r1, r2, r3, r4, r5, r6]
            num_possible_products_list = [1, 1, 2, 1, 2, 1]

            i = 0
            for i in range(6):
                if alkene_regiochemistry_list[i] == True:
                    j = 0
                    for j in range(num_possible_products_list[i]):
                        smile_products.append(
                            Chem.MolToSmiles(H_halideaddition_products_list[i][j][0]))
                        j += 1
                    del (smile_products[0])
                    smile1 = smile_products[0]

                    i = 0
                    break

                i += 1

        k = 0
        # for k in range(len(smile_products)):
        #     draw(smile_products[k])
        #     k += 1

        return smile_products

    @staticmethod
    def catalytic_hydrogenation(smile):
        smile_products = [smile]
        # draw(smile)

        smile1 = smile_products[0]

        while FunctionalGroups.alkene(smile1) == True:
            rxn1 = AllChem.ReactionFromSmarts('[C:1]=[C:2].[H:3][H:4]>>[C:1](-[H:3])[C:2](-[H:4])',
                                              useSmiles=True)
            ps1 = rxn1.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H][H]")))
            smile_products.append(Chem.MolToSmiles(ps1[0][0]))
            del (smile_products[0])
            smile1 = smile_products[0]

        k = 0
        # for k in range(len(smile_products)):
        #     draw(smile_products[k])
        #     k += 1

        return smile_products

    @staticmethod
    def epoxidation(smile):
        smile_products = [smile]
        # draw(smile)

        smile1 = smile_products[0]

        while FunctionalGroups.alkene(smile1) == True:
            rxn1 = AllChem.ReactionFromSmarts('[C:1]=[C:2].[O:3]>>[C:1]1[C:2][O:3]1',
                                              useSmiles=True)
            ps1 = rxn1.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("O")))
            smile_products.append(Chem.MolToSmiles(ps1[0][0]))
            del (smile_products[0])
            smile1 = smile_products[0]

        k = 0
        # for k in range(len(smile_products)):
        #     draw(smile_products[k])
        #     k += 1

        return smile_products

    @staticmethod
    def hydroxylation(smile):
        smile_products = [smile]
        # draw(smile)

        smile1 = smile_products[0]

        while FunctionalGroups.epoxide(smile1) == True:
            rxn1 = AllChem.ReactionFromSmarts('[C:1]1[C:2][O:3]1.[O:4]>>[C:1]([O:3])[C:2]([O:4])',
                                              useSmiles=True)
            ps1 = rxn1.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("O")))
            smile_products.append(Chem.MolToSmiles(ps1[0][0]))
            del (smile_products[0])
            smile1 = smile_products[0]

        k = 0
        # for k in range(len(smile_products)):
        #     draw(smile_products[k])
        #     k += 1

        return smile_products

    @staticmethod
    def dichlorocyclopropane_formation(smile):
        smile_products = [smile]
        # draw(smile)

        smile1 = smile_products[0]

        while FunctionalGroups.alkene(smile1) == True:
            rxn1 = AllChem.ReactionFromSmarts(
                '[C:1]=[C:2].[Cl:3][C:4][Cl:5]>>[Cl:3][C:4]1([Cl:5])[C:1][C:2]1', useSmiles=True)
            ps1 = rxn1.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("ClCCl")))
            smile_products.append(Chem.MolToSmiles(ps1[0][0]))
            del (smile_products[0])
            smile1 = smile_products[0]

        k = 0
        # for k in range(len(smile_products)):
        #     draw(smile_products[k])
        #     k += 1

        return smile_products

    @staticmethod
    def dibromocyclopropane_formation(smile):
        smile_products = [smile]
        # draw(smile)

        smile1 = smile_products[0]

        while FunctionalGroups.alkene(smile1) == True:
            rxn1 = AllChem.ReactionFromSmarts(
                '[C:1]=[C:2].[Br:3][C:4][Br:5]>>[Br:3][C:4]1([Br:5])[C:1][C:2]1', useSmiles=True)
            ps1 = rxn1.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("BrCBr")))
            smile_products.append(Chem.MolToSmiles(ps1[0][0]))
            del (smile_products[0])
            smile1 = smile_products[0]

        k = 0
        # for k in range(len(smile_products)):
        #     draw(smile_products[k])
        #     k += 1
        return smile_products

    @staticmethod
    def diiodocyclopropane_formation(smile):
        smile_products = [smile]
        # draw(smile)

        smile1 = smile_products[0]

        while FunctionalGroups.alkene(smile1) == True:
            rxn1 = AllChem.ReactionFromSmarts(
                '[C:1]=[C:2].[I:3][C:4][I:5]>>[I:3][C:4]1([I:5])[C:1][C:2]1', useSmiles=True)
            ps1 = rxn1.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("ICI")))
            smile_products.append(Chem.MolToSmiles(ps1[0][0]))
            del (smile_products[0])
            smile1 = smile_products[0]

        k = 0
        # for k in range(len(smile_products)):
        #     draw(smile_products[k])
        #     k += 1
        return smile_products

    @staticmethod
    def Simmons_Smith_reaction(smile):
        smile_products = [smile]
        # draw(smile)

        smile1 = smile_products[0]

        while FunctionalGroups.alkene(smile1) == True:
            rxn1 = AllChem.ReactionFromSmarts('[C:1]=[C:2].[C:3]>>[C:1]1[C:2][C:3]1',
                                              useSmiles=True)
            ps1 = rxn1.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("C")))
            smile_products.append(Chem.MolToSmiles(ps1[0][0]))
            del (smile_products[0])
            smile1 = smile_products[0]

        k = 0
        # for k in range(len(smile_products)):
        #     draw(smile_products[k])
        #     k += 1
        return smile_products

    @staticmethod
    def alkene_ozonation(smile):
        # draw(smile)
        smile_products = [smile]

        i = 0
        while i in range(len(smile_products)):
            if FunctionalGroups.ring_alkene(smile_products[i]) == True:
                # if alkene is in ring, ensure that the cleavage results in one molecule, not two
                m = Chem.MolFromSmiles(smile_products[i])

                bis = m.GetSubstructMatches(Chem.MolFromSmarts('[C;R]=[C;R]'))
                bs = [m.GetBondBetweenAtoms(x, y).GetIdx() for x, y in bis]
                bs = []
                labels = []
                for bi in bis:
                    b = m.GetBondBetweenAtoms(bi[0], bi[1])
                    if b.GetBeginAtomIdx() == bi[0]:
                        labels.append((100, 100))
                    else:
                        labels.append((100, 100))
                    bs.append(b.GetIdx())

                nm = Chem.FragmentOnBonds(m, bs, dummyLabels=labels)
                smile_products.append(Chem.MolToSmiles(nm, True))
                del (smile_products[i])

                while FunctionalGroups.cleaved_alkene(smile_products[i]) == True:
                    rxn1 = AllChem.ReactionFromSmarts(
                        '[C:1]=[100*:2].[O:3]>>[C:1](=[O:3]).[100*:2]', useSmiles=True)
                    ps1 = rxn1.RunReactants(
                        (Chem.MolFromSmiles(smile_products[i]), Chem.MolFromSmiles("O")))
                    smile_products.append(Chem.MolToSmiles(ps1[0][0]))
                    del (smile_products[i])
                    i = i

            elif FunctionalGroups.alkene(smile_products[i]) == True:
                # alkene cleavage should result in two molecules
                rxn1 = AllChem.ReactionFromSmarts('[C:1]=[C:2].[O:3][O:4]>>[C:1]=[O:3].[C:2]=[O:4]',
                                                  useSmiles=True)
                ps1 = rxn1.RunReactants(
                    (Chem.MolFromSmiles(smile_products[i]), Chem.MolFromSmiles("[O][O]")))
                smile_products.append(Chem.MolToSmiles(ps1[0][0]))
                smile_products.append(Chem.MolToSmiles(ps1[0][1]))
                del (smile_products[i])
                i = 0

            else:
                i += 1

        k = 0
        # for k in range(len(smile_products)):
        #     draw(smile_products[k])
        #     k += 1
        return smile_products

    # @staticmethod
    # def KMnO4_alkene_cleavage(smile):
    #     # draw(smile)
    #     smile_products = [smile]
    #
    #     i = 0
    #     while i in range(len(smile_products)):
    #         if FunctionalGroups.ring_alkene(smile_products[i]) == True:
    #             m = Chem.MolFromSmiles(smile_products[i])
    #
    #             bis = m.GetSubstructMatches(Chem.MolFromSmarts('[C;R]=[C;R]'))
    #             bs = [m.GetBondBetweenAtoms(x, y).GetIdx() for x, y in bis]
    #             bs = []
    #             labels = []
    #             for bi in bis:
    #                 b = m.GetBondBetweenAtoms(bi[0], bi[1])
    #                 if b.GetBeginAtomIdx() == bi[0]:
    #                     labels.append((100, 100))
    #                 else:
    #                     labels.append((100, 100))
    #                 bs.append(b.GetIdx())
    #
    #             nm = Chem.FragmentOnBonds(m, bs, dummyLabels=labels)
    #             smile_products.append(Chem.MolToSmiles(nm, True))
    #             print(Chem.MolToSmiles(nm, True))
    #             del (smile_products[i])
    #
    #             while FunctionalGroups.cleaved_alkene(smile_products[i]) == True:
    #                 # One hydrogen bonded to the carbon
    #                 rxn1 = AllChem.ReactionFromSmarts(
    #                     '[C:0][CH:1]=[100*:2].[O:3]=[O:4]>>[C:0][C:1](=[O:3])-[O:4].[100*:2]',
    #                     useSmiles=True)
    #                 ps1 = rxn1.RunReactants(
    #                     (Chem.MolFromSmiles(smile_products[i]), Chem.MolFromSmiles("O=O")))
    #
    #                 # No hydrogens bonded to the carbon
    #                 rxn2 = AllChem.ReactionFromSmarts('[C:1]=[100*:2].[O:3]>>[C:1]=[O:3].[100*:2]',
    #                                                   useSmiles=True)
    #                 ps2 = rxn2.RunReactants(
    #                     (Chem.MolFromSmiles(smile_products[i]), Chem.MolFromSmiles("O")))
    #
    #                 m1 = Chem.MolFromSmiles(smile_products[i])
    #                 patt1 = Chem.MolFromSmarts("[CH]=[100*]")
    #                 r1 = m1.HasSubstructMatch(patt1)
    #
    #                 m2 = Chem.MolFromSmiles(smile_products[i])
    #                 patt2 = Chem.MolFromSmarts("[CH0](=[100*])")
    #                 r2 = m2.HasSubstructMatch(patt2)
    #
    #                 if r1 == True:
    #                     smile_products.append(Chem.MolToSmiles(ps1[0][0]))
    #                     del (smile_products[i])
    #                     i = 0
    #                 elif r2 == True:
    #                     smile_products.append(Chem.MolToSmiles(ps2[0][0]))
    #                     del (smile_products[i])
    #                     i = 0
    #
    #         elif FunctionalGroups.nonring_alkene(smile_products[i]) == True:
    #
    #             # Two hydrogens bonded to the carbon
    #             rxn1 = AllChem.ReactionFromSmarts(
    #                 '[CH2:1]=[CH2:2].([O:3]=[O:4].[O:5]=[O:6])>>[O:3]=[C:1]=[O:4].[O:5]=[C:2]=[O:6]',
    #                 useSmiles=True)
    #             ps1 = rxn1.RunReactants(
    #                 (Chem.MolFromSmiles(smile_products[i]), Chem.MolFromSmiles("O=O.O=O")))
    #
    #             # One hydrogen bonded to the carbon
    #             rxn2 = AllChem.ReactionFromSmarts(
    #                 '[CH2:1]=[CH:2].([O:3]=[O:4].[O:5]=[O:6])>>[O:3]=[C:1]=[O:4].[C:2]([O:5])=[O:6]',
    #                 useSmiles=True)
    #             ps2 = rxn2.RunReactants(
    #                 (Chem.MolFromSmiles(smile_products[i]), Chem.MolFromSmiles("O=O.O=O")))
    #
    #             # No hydrogens bonded to the carbon
    #             rxn3 = AllChem.ReactionFromSmarts(
    #                 '[CH2:1]=[CH0:2].([O:3]=[O:4].[O:5])>>[O:3]=[C:1]=[O:4].[CH0:2]=[O:5]',
    #                 useSmiles=True)
    #             ps3 = rxn3.RunReactants(
    #                 (Chem.MolFromSmiles(smile_products[i]), Chem.MolFromSmiles("O=O.O")))
    #
    #             rxn4 = AllChem.ReactionFromSmarts(
    #                 '[CH:1]=[CH:2].([O:3]=[O:4].[O:5]=[O:6])>>[C:1]([O:3])=[O:4].[C:2]([O:5])=[O:6]',
    #                 useSmiles=True)
    #             ps4 = rxn4.RunReactants(
    #                 (Chem.MolFromSmiles(smile_products[i]), Chem.MolFromSmiles("O=O.O=O")))
    #
    #             rxn5 = AllChem.ReactionFromSmarts(
    #                 '[CH:1]=[CH0:2].([O:3]=[O:4].[O:5])>>[C:1]([O:3])=[O:4].[CH0:2]=[O:5]',
    #                 useSmiles=True)
    #             ps5 = rxn5.RunReactants(
    #                 (Chem.MolFromSmiles(smile_products[i]), Chem.MolFromSmiles("O=O.O")))
    #
    #             rxn6 = AllChem.ReactionFromSmarts(
    #                 '[CH0:1]=[CH0:2].[O:3]=[O:4]>>[CH0:1]=[O:3].[CH0:2]=[O:4]', useSmiles=True)
    #             ps6 = rxn6.RunReactants(
    #                 (Chem.MolFromSmiles(smile_products[i]), Chem.MolFromSmiles("O=O")))
    #
    #             m1 = Chem.MolFromSmiles(smile_products[i])
    #             patt1 = Chem.MolFromSmarts("[CH2]=[CH2]")
    #             r1 = m1.HasSubstructMatch(patt1)
    #
    #             m2 = Chem.MolFromSmiles(smile_products[i])
    #             patt2 = Chem.MolFromSmarts("[CH2]=[CH]")
    #             r2 = m2.HasSubstructMatch(patt2)
    #
    #             m3 = Chem.MolFromSmiles(smile_products[i])
    #             patt3 = Chem.MolFromSmarts("[CH2]=[CH0]")
    #             r3 = m3.HasSubstructMatch(patt3)
    #
    #             m4 = Chem.MolFromSmiles(smile_products[i])
    #             patt4 = Chem.MolFromSmarts("[CH]=[CH]")
    #             r4 = m4.HasSubstructMatch(patt4)
    #
    #             m5 = Chem.MolFromSmiles(smile_products[i])
    #             patt5 = Chem.MolFromSmarts("[CH]=[CH0]")
    #             r5 = m5.HasSubstructMatch(patt5)
    #
    #             m6 = Chem.MolFromSmiles(smile_products[i])
    #             patt6 = Chem.MolFromSmarts("[CH0]=[CH0]")
    #             r6 = m6.HasSubstructMatch(patt6)
    #
    #             list_of_possible_configurations = [r1, r2, r3, r4, r5, r6]
    #             list_of_possible_reactions = [ps1, ps2, ps3, ps4, ps5, ps6]
    #
    #             j = 0
    #             for j in range(6):
    #                 if list_of_possible_configurations[j] == True:
    #                     smile_products.append(Chem.MolToSmiles(list_of_possible_reactions[j][0][0]))
    #                     smile_products.append(Chem.MolToSmiles(list_of_possible_reactions[j][0][1]))
    #
    #                     del (smile_products[i])
    #                     i = i
    #
    #                 else:
    #                     j += 1
    #
    #         else:
    #             i += 1
    #
    #     k = 0
    #     # for k in range(len(smile_products)):
    #     #     draw(smile_products[k])
    #     #     k += 1
    #
    #     return smile_products

    @staticmethod
    def oxidative_cleavage_of_1_2_diols(smile):
        # draw(smile)
        smile_products = [smile]

        i = 0
        while i in range(len(smile_products)):
            if FunctionalGroups.ring_1_2_diols(smile_products[i]) == True:

                m = Chem.MolFromSmiles(smile_products[i])
                bis = m.GetSubstructMatches(Chem.MolFromSmarts('[C;R]([O])[C;R]([O])'))

                bs = [m.GetBondBetweenAtoms(a, c).GetIdx() for a, b, c, d in bis]

                bs = []
                labels = []
                for bi in bis:
                    b = m.GetBondBetweenAtoms(bi[0], bi[2])
                    if b.GetBeginAtomIdx() == bi[0]:
                        labels.append((100, 100))
                    else:
                        labels.append((100, 100))
                    bs.append(b.GetIdx())

                nm = Chem.FragmentOnBonds(m, bs, dummyLabels=labels)
                smile_products.append(Chem.MolToSmiles(nm, True))
                del (smile_products[i])

                while FunctionalGroups.cleaved_ring_vicinal_diols(smile_products[i]) == True:
                    rxn1 = AllChem.ReactionFromSmarts(
                        '[C:1]([100*:2])[O:3].[O:4]>>[C:1]=[O:4].[100*:2].[O:3]', useSmiles=True)
                    ps1 = rxn1.RunReactants(
                        (Chem.MolFromSmiles(smile_products[i]), Chem.MolFromSmiles("O")))
                    smile_products.append(Chem.MolToSmiles(ps1[0][0]))
                    del (smile_products[i])
                    i = i

            elif FunctionalGroups.non_ring_1_2_diols(smile_products[i]) == True:
                rxn1 = AllChem.ReactionFromSmarts(
                    '[C:1]([O:2])-[C:3]([O:4]).[O:5][O:6]>>[C:1]=[O:5].[C:3]=[O:6].[O:2][O:4]',
                    useSmiles=True)
                ps1 = rxn1.RunReactants(
                    (Chem.MolFromSmiles(smile_products[i]), Chem.MolFromSmiles("[O][O]")))
                smile_products.append(Chem.MolToSmiles(ps1[0][0]))
                smile_products.append(Chem.MolToSmiles(ps1[0][1]))
                del (smile_products[i])
                i = 0

            else:
                i += 1

        k = 0
        # for k in range(len(smile_products)):
        #     draw(smile_products[k])
        #     k += 1
        return smile_products

    @staticmethod
    def partial_alkyne_hydrogenation(smile):
        smile_products = [smile]
        # draw(smile)

        smile1 = smile_products[0]

        while FunctionalGroups.alkyne(smile1) == True:
            rxn1 = AllChem.ReactionFromSmarts('[C:1]#[C:2].[H:3][H:4]>>[C:1](-[H:3])=[C:2](-[H:4])',
                                              useSmiles=True)
            ps1 = rxn1.RunReactants((Chem.MolFromSmiles(smile1), Chem.MolFromSmiles("[H][H]")))
            smile_products.append(Chem.MolToSmiles(ps1[0][0]))
            del (smile_products[0])
            smile1 = smile_products[0]

        # k = 0
        # for k in range(len(smile_products)):
        #     draw(smile_products[k])
        #     k += 1

        return smile_products



    @staticmethod
    def getAllReactions():
        method_list = [eval("AlkeneReactions." + method) for method in dir(AlkeneReactions) if
                       method.startswith('_') is False and method.startswith("getAll") is False]
        return method_list
