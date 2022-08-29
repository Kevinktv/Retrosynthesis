from rdkit import Chem


class FunctionalGroups:


    # UPDATED METHOD: This may possibly much much simplify the problem of finding specific functional groups in a given molecule

    @staticmethod
    def alkene(smile):
        m = Chem.MolFromSmiles(smile)
        patt = Chem.MolFromSmarts('C=C')
        FG_found = m.HasSubstructMatch(patt)

        if FG_found == True:
            return True
        else:
            return False

    @staticmethod
    def alkyne(smile):
        m = Chem.MolFromSmiles(smile)
        patt = Chem.MolFromSmarts('C#C')
        FG_found = m.HasSubstructMatch(pat)

        if FG_found == True:
            return True
        else:
            return False

    @staticmethod
    def benzene(smile):
        m = Chem.MolFromSmiles(smile)
        patt = Chem.MolFromSmarts("c1ccccc1")
        FG_found = m.HasSubstructMatch(patt)

        if FG_found == True:
            return True
        else:
            return False

    @staticmethod
    def amine(smile):
        m = Chem.MolFromSmiles(smile)
        patt1 = Chem.MolFromSmarts("NC")
        patt2 = Chem.MolFromSmarts("Nc")
        FG_found1 = m.HasSubstructMatch(patt1)
        FG_found2 = m.HasSubstructMatch(patt2)

        if FG_found1 == True or FG_found2 == True:
            return True
        else:
            return False

    @staticmethod
    def alcohol(smile):
        m = Chem.MolFromSmiles(smile)
        patt1 = Chem.MolFromSmiles("O[H]")
        FG_found1 = m.GetSubstructMatches(patt1)

        patt2 = Chem.MolFromSmiles("[H]OC=O")
        FG_found2 = m.GetSubstructMatches(patt2)

        if len(FG_found1) > len(FG_found2):
            return True
        else:
            return False

    @staticmethod
    def ether(smile):
        m = Chem.MolFromSmiles(smile)
        patt = Chem.MolFromSmiles("COC")

        FG_found = m.HasSubstructMatch(patt)

        if FG_found == True:
            return True
        else:
            return False

    @staticmethod
    def aldehyde(smile):
        m = Chem.MolFromSmiles(smile)
        patt = Chem.MolFromSmiles("O=CC")

        FG_found = m.HasSubstructMatch(patt)

        if FG_found == True:
            return True
        else:
            return False

    @staticmethod
    def ketone(smile):
        m = Chem.MolFromSmiles(smile)
        patt = Chem.MolFromSmiles("CC(C)=O")

        FG_found = m.HasSubstructMatch(patt)

        if FG_found == True:
            return True
        else:
            return False

    @staticmethod
    def alkylhalide(smile):
        m = Chem.MolFromSmiles(smile)
        patt1 = Chem.MolFromSmiles("Br")
        patt2 = Chem.MolFromSmiles("Cl")
        patt3 = Chem.MolFromSmiles("F")
        patt4 = Chem.MolFromSmiles("I")

        FG_found1 = m.HasSubstructMatch(patt1)
        FG_found2 = m.HasSubstructMatch(patt2)
        FG_found3 = m.HasSubstructMatch(patt3)
        FG_found4 = m.HasSubstructMatch(patt4)

        if FG_found1 == True or FG_found2 == True or FG_found3 == True or FG_found4 == True:
            return True
        else:
            return False

    @staticmethod
    def thiol(smile):
        m = Chem.MolFromSmiles(smile)
        patt = Chem.MolFromSmiles("CS[H]")
        FG_found = m.HasSubstructMatch(patt)

        if FG_found == True:
            return True
        else:
            return False

    @staticmethod
    def carboxylicacid(smile):
        m = Chem.MolFromSmiles(smile)
        patt = Chem.MolFromSmiles("CC(O)=O")
        FG_found = m.HasSubstructMatch(patt)

        if FG_found == True:
            return True
        else:
            return False

    @staticmethod
    def ester(smile):
        m = Chem.MolFromSmiles(smile)
        patt = Chem.MolFromSmiles("CC(OC)=O")
        FG_found = m.HasSubstructMatch(patt)

        if FG_found == True:
            return True
        else:
            return False

    @staticmethod
    def amide(smile):
        m = Chem.MolFromSmiles(smile)
        patt = Chem.MolFromSmiles("C(N)=O")
        FG_found = m.HasSubstructMatch(patt)

        if FG_found == True:
            return True
        else:
            return False
