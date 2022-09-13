"""
This file contains various attempts to make a retrosynthesis AI.
Yay
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
                                                AlkeneReactions.getAllReactions(), 2)

    closestSimilarity = -1
    indexOfClosest = -1
    print(products)

    for index in products:
        for product in products[index]:
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
    minimaxhelper(startSmile, endingSmile, possibleReactions, bestReactionSet, [], depth)
    return bestReactionSet[0]

def minimaxhelper(startSmile, endingSmile, possibleReactions,
                  bestReactionSet, currentReacSet, depth):
    if startSmile == endingSmile or depth == 0:
        return evaluate(startSmile, endingSmile)

    maxEval = -999 # Initially set eval of that posn to -999
    for reaction in possibleReactions:
        currentReacSet.append(reaction)
        # Carry out the reaction and then evaluate the new smile with ending smile
        newSmileList = reaction(startSmile)

        for newSmile in newSmileList:
            val = minimaxhelper(newSmile, endingSmile, possibleReactions, bestReactionSet, currentReacSet, depth - 1)
            if val > maxEval:
                maxEval = val
                bestReactionSet[0] = currentReacSet[:]

        currentReacSet.pop()

    return maxEval

# def minimaxhelper(startSmile, endingSmile, possibleReactions,
#                   bestReactionSet, currentReacSet, bestVal, depth):
#     if startSmile == endingSmile or depth == 0:
#         val = evaluate(startSmile, endingSmile)
#         if val > bestVal[0]:
#             bestReactionSet[0] = currentReacSet[:]
#             bestVal[0] = val
#         return
#
#     print("1 " ,currentReacSet)
#     for reaction in possibleReactions:
#         currentReacSet.append(reaction)
#         # Carry out the reaction and then evaluate the new smile with ending smile
#         newSmileList = reaction(startSmile)
#         print("2", newSmileList)
#         for newSmile in newSmileList:
#             minimaxhelper(newSmile, endingSmile, possibleReactions, bestReactionSet,
#                           currentReacSet, bestVal, depth - 1)
#
#         currentReacSet.pop()



# ###########################
# from copy import deepcopy
# from montecarlo.node import Node
# from montecarlo.montecarlo import MonteCarlo
#
#
# class GameNode:
#     def __init__(self, currsmile: str, requiredsmile: str):
#         self.smile = currsmile
#         self.endsmile = requiredsmile
#         self.children = []
#
#     def getPossibleMoves(self):
#         # Return all possible moves/reactions from current point
#         # return AlkeneReactions.getAllReactions()
#         return [Reactions.florination, Reactions.chlorination]
#
#     def move(self, reaction):
#         # Makes a reaction and updates the state. Reaction is of type function
#         self.smile = reaction(self.smile)
#
#     def addChild(self, child):
#         # Adds child gamenode to children list
#         self.children.append(child)
#
#     def evaluate(self):
#         ms = [Chem.MolFromSmiles(self.smile), Chem.MolFromSmiles(self.endsmile)]
#
#         # for example
#
#         fps = [Chem.RDKFingerprint(x) for x in ms]
#         return DataStructs.FingerprintSimilarity(fps[0], fps[1])
#
#
# def child_finder(self, node):
#     for move in node.state.getPossibleMoves():
#         child = Node(deepcopy(node.state))  # Make copy of parent node
#         child.state.move(move)  # or however your library works
#         node.add_child(child)
#
#
# def node_evaluator(self, node):
#     if node.state.won():
#         return 1
#     elif node.state.lost():
#         return -1
#
#
# game = GameNode("c1ccccc1", "c1c")
#
# montecarlo = MonteCarlo(Node(game))
# montecarlo.child_finder = child_finder
# montecarlo.node_evaluator = node_evaluator

########################
import numpy as np
from collections import defaultdict


class MonteCarloTreeSearchNode():
    def __init__(self, state, parent=None, parent_action=None):
        self.state = state
        self.parent = parent
        self.parent_action = parent_action
        self.children = []
        self._number_of_visits = 0
        self._results = defaultdict(int)
        self._results[1] = 0
        self._results[-1] = 0
        self._untried_actions = None
        self._untried_actions = self.untried_actions()

        # self.final_product = final_product
        return

    def untried_actions(self):
        self._untried_actions = self.state.get_legal_actions()
        return self._untried_actions

    def q(self):
        wins = self._results[1]
        loses = self._results[-1]
        return wins - loses

    def n(self):
        return self._number_of_visits

    def expand(self):
        action = self._untried_actions.pop()
        next_state = self.state.move(action)
        child_node = MonteCarloTreeSearchNode(
            next_state, parent=self, parent_action=action)

        self.children.append(child_node)
        return child_node

    def is_terminal_node(self):
        return self.state.is_game_over()

    def rollout(self):
        current_rollout_state = self.state

        while not current_rollout_state.is_game_over():
            possible_moves = current_rollout_state.get_legal_actions()

            action = self.rollout_policy(possible_moves)
            current_rollout_state = current_rollout_state.move(action)
        return current_rollout_state.game_result()

    def backpropagate(self, result):
        self._number_of_visits += 1.
        self._results[result] += 1.
        if self.parent:
            self.parent.backpropagate(result)

    def is_fully_expanded(self):
        return len(self._untried_actions) == 0

    def best_child(self, c_param=0.1):

        choices_weights = [(c.q() / c.n()) + c_param * np.sqrt((2 * np.log(self.n()) / c.n())) for c
                           in self.children]
        return self.children[np.argmax(choices_weights)]

    def rollout_policy(self, possible_moves):

        return possible_moves[np.random.randint(len(possible_moves))]

    def _tree_policy(self):

        current_node = self
        while not current_node.is_terminal_node():

            if not current_node.is_fully_expanded():
                return current_node.expand()
            else:
                current_node = current_node.best_child()
        return current_node

    def best_action(self):
        simulation_no = 100

        for i in range(simulation_no):
            v = self._tree_policy()
            reward = v.rollout()
            v.backpropagate(reward)

        return self.best_child(c_param=0.)

class State:
    def __init__(self, start, end):
        self.start = start
        self.end = end


    def get_legal_actions(self):
        '''
        Modify according to your game or
        needs. Constructs a list of all
        possible actions from current state.
        Returns a list.
        '''
        return [Reactions.chlorination, Reactions.florination]

    def evaluate(self, smile1, smile2):
        """
        Evaluate two smiles for their similarity and returns how similar they are.
        """
        # Evaluates and returns how similar smile1 and smile2 are
        ms = [Chem.MolFromSmiles(smile1), Chem.MolFromSmiles(smile2)]

        # for example

        fps = [Chem.RDKFingerprint(x) for x in ms]
        return DataStructs.FingerprintSimilarity(fps[0], fps[1])

    def is_game_over(self):
        '''
        It is the game over condition
        and depends on your game. Returns
        true or false
        '''
        evalu = self.evaluate(self.start, self.end)
        if evalu >= 0.7:
            return True

        if evalu <= 0.3:
            return False

    def game_result(self):
        '''
        Returns 1 or 0 or -1 depending
        on your state corresponding to win,
        tie or a loss.
        '''
        return self.evaluate(self.start, self.end)

    def move(self, action):
        '''
        '''
        self.start = action(self.start)


def main():
    initial_state = State("CC=C", "CCCC")
    finalprod = "Clc1ccccc1"
    root = MonteCarloTreeSearchNode(state=initial_state)
    selected_node = root.best_action()
    return
