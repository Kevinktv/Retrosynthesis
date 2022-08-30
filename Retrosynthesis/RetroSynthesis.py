"""
This file contains various attempts to make a retrosynthesis AI
"""
from rdkit import Chem
from rdkit import DataStructs

#######################################
from Retrosynthesis.ChemReac.reactions import AlkeneReactions, Reactions
#######################################


def allPossibleReactions(reactions: list, depth: int):
    """
    This function generates a list containing every combination of the reaction in reactions
    till depth.

    :param reactions:  List containing various reactions
    :param depth: The no of reactions in each combination
    :return: Returns a list containing all reaction combinations
    """
    res = []

    allPossibleReactionsHelper(reactions, res, [], depth)
    return res


def allPossibleReactionsHelper(reactions, res, path, depth):
    """
    Helper for above method
    """
    if len(path) == depth:
        res.append(path[:])
        return

    for reaction in reactions:
        path.append(reaction)
        allPossibleReactionsHelper(reactions, res, path, depth)
        path.pop()


def carryOutReactions(smile, reactions, depth):
    """
    Given a list containing various reactions, this method evaluates all combinations of
    those reactions till depth and then runs them on smile. It returns a tuple where index 0
    are the reactions and index 1 is a dictionary that maps the reactions index to the products
    formed by that set of reactions.

    :param smile: Starting smile
    :param reactions: A list of reactions
    :param depth: The depth of reactions formed for combinations.
    :return:
    """
    reactionLists = allPossibleReactions(reactions, depth)  # Getting all possible permutations
    res = {}
    index = 0
    for reactionList in reactionLists:  # Every individual permutation of reactions
        products = []
        for reaction in reactionList:  # Every reaction in that particular permutation
            if products == []:
                products = reaction(smile)
            else:
                temp = []
                for product in products:
                    temp.extend(reaction(product))

                products = temp

        res[index] = products
        index += 1

    return (reactionLists, res)


def evaluate(smile1, smile2):
    """
    Evaluate two smiles for their similarity and returns how similar they are.
    """
    # Evaluates and returns how similar smile1 and smile2 are
    ms = [Chem.MolFromSmiles(smile1), Chem.MolFromSmiles(smile2)]

    # for example

    fps = [Chem.RDKFingerprint(x) for x in ms]
    return DataStructs.FingerprintSimilarity(fps[0], fps[1])


def naiveRetroSynthesis(start, end):
    """
    In the naive retrosynthesis, we are simply doing every possible reaction combinations and then
    running evaluate on each of them to see how close it is to "end". We then return the reactions
    that cause it.

    Currently I am using only chlorination and florination

    start: The starting smile
    end: The finishing smile we desire
    :return:
    """
    reactionOrder, products = carryOutReactions(start,
                                                [Reactions.chlorination, Reactions.florination], 2)

    closestSimilarity = -1
    indexOfClosest = -1

    for index in products:
        for product in products[index]:
            print("######")
            print(product)
            print(end)
            similarity = evaluate(end, product)
            if similarity > closestSimilarity:
                indexOfClosest = index
                closestSimilarity = similarity

    return reactionOrder[indexOfClosest]
################
# Below, I am making a minimax based approach of trying to choose the best move after
# backtracking rather than creating all possible combinations.

def minimax(startSmile, endingSmile, possibleReactions, depth):

    bestReactionSet = [[]]
    minimaxhelper(startSmile, endingSmile, possibleReactions, bestReactionSet, [], [-999], depth)
    return bestReactionSet[0]

# def minimaxhelper(startSmile, endingSmile, possibleReactions,
#                   bestReactionSet, currentReacSet, depth):
#     if startSmile == endingSmile or depth == 0:
#         return evaluate(startSmile, endingSmile)
#
#     maxEval = -999 # Initially set eval of that posn to -999
#     for reaction in possibleReactions:
#         currentReacSet.append(reaction)
#         # Carry out the reaction and then evaluate the new smile with ending smile
#         newSmile = reaction(startSmile)
#         val = minimaxhelper(newSmile, endingSmile, possibleReactions, depth - 1)
#         if val > maxEval:
#             maxEval = val
#             bestReactionSet = currentReacSet[:]
#
#         currentReacSet.pop()
#
#     return maxEval

def minimaxhelper(startSmile, endingSmile, possibleReactions,
                  bestReactionSet, currentReacSet, bestVal, depth):
    if startSmile == endingSmile or depth == 0:
        val = evaluate(startSmile, endingSmile)
        if val > bestVal[0]:
            bestReactionSet[0] = currentReacSet[:]
            bestVal[0] = val
        return

    for reaction in possibleReactions:
        currentReacSet.append(reaction)
        # Carry out the reaction and then evaluate the new smile with ending smile
        newSmileList = reaction(startSmile)
        for newSmile in newSmileList:
            minimaxhelper(newSmile, endingSmile, possibleReactions, bestReactionSet,
                          currentReacSet, bestVal, depth - 1)

        currentReacSet.pop()
