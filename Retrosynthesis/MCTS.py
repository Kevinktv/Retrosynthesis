from __future__ import print_function

import random
import time
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit import DataStructs

# Needed for MCTS
from collections import defaultdict
import numpy as np

#######
from Retrosynthesis.ChemReac.functional_groups import FunctionalGroups
from Retrosynthesis.ChemReac.reactions import AlkeneReactions

# -*- coding: utf-8 -*-
"""Copy of MCTS Attempt

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1uAu6QAhI-PKTFA-dP1z7-MDlF9aL4hzH

# Monte Carlo Tree Search

## The connect4 game class
"""

def draw(smile):
    molecule = Chem.MolFromSmiles(smile)
    fig = Draw.MolToMPL(molecule, size=(200, 200))


class Game:

    def __init__(self, initial_SMILE, final_SMILE):
        self.board = initial_SMILE
        self.end = final_SMILE
        self.turn = 1
        self.history = []
        self.ctr = 0

    # check whether or not the game is over
    def check_win(self):
        # Evaluates and returns how similar smile1 and smile2 are
        ms = [Chem.MolFromSmiles(self.board), Chem.MolFromSmiles(self.end)]

        fps = [Chem.RDKFingerprint(x) for x in ms]
        return DataStructs.FingerprintSimilarity(fps[0], fps[1]) == 1 or self.ctr >= 10

    # make a move in column x
    def make_move(self, reaction):
        self.history.append(self.board)
        self.board = reaction(self.board)[0]
        self.ctr += 1
        return self.check_win()

    # unmake the last move in column x
    def unmake_move(self):
        self.board = self.history.pop()
        self.ctr -= 1

    # return a list of available moves
    def moves(self):
        return AlkeneReactions.getAllReactions()

    # print the board
    def show(self):
        print("THE BOARD IS: ", self.board, self.history)

    def check_win2(self):
        ms = [Chem.MolFromSmiles(self.board), Chem.MolFromSmiles(self.end)]

        fps = [Chem.RDKFingerprint(x) for x in ms]
        return DataStructs.FingerprintSimilarity(fps[0], fps[1]) == 1

    def evaluate(self):
        ms = [Chem.MolFromSmiles(self.board), Chem.MolFromSmiles(self.end)]

        fps = [Chem.RDKFingerprint(x) for x in ms]
        return DataStructs.FingerprintSimilarity(fps[0], fps[1])


"""## Alpha-beta search"""


class AlphaBeta:
    def __init__(self, initial_SMILE, final_SMILE, depth=8):
        self.game = Game(initial_SMILE, final_SMILE)
        self.depth = depth

    # returns a move based on an alpha-beta search
    def act(self):
        depth = 9
        move = self.search()
        return move

    # update the internal board state for the class
    def feed(self, move):
        board = self.game.board
        end = self.game.end
        newboard = move(board)[0]
        self.game = Game(newboard, end)
        print(self.game.board)

    # the root node of an alpha-beta search
    def search(self):
        print("AlphaBeta searching")
        moves = self.game.moves()

        # a list to store the values associated with each move
        scores = []
        alpha = -10
        for move in moves:
            res = self.game.make_move(move)
            # if the move wins the game, play it immediately

            if res:
                self.game.unmake_move()
                return move
            # val = -self.alpha_beta(-10, -alpha, self.depth - 1)
            val = self.alpha_beta(self.depth-1)
            self.game.unmake_move()
            scores.append((val, move))


        # the algorithm randomises between moves that have the same value
        random.shuffle(scores)
        scores.sort(key=lambda x: -x[0])
        print("AlphaBeta score: " + str(scores[0][0]))
        return scores[0][1]

    # def alpha_beta(self, alpha, beta, depth):
    #     print("broooo", depth)
    #     if depth == 0: return 0
    #     moves = self.game.moves()
    #     # check whether or not the game is drawn
    #     if self.game.check_win(): return 0
    #     for move in moves:
    #         res = self.game.make_move(move)
    #         # if the move wins the game, return a winning score
    #         if res:
    #             self.game.unmake_move()
    #             return 1 + 0.01 * depth
    #         val = -self.alpha_beta(-beta, -alpha, depth - 1)
    #         self.game.unmake_move()
    #
    #         # check for alpha node
    #         if val >= alpha:
    #             alpha = val
    #
    #         # check for beta cut
    #         if val >= beta:
    #             return val
    #     return alpha
    def alpha_beta(self, depth):
        if self.game.check_win2() or depth == 0:
            return self.game.evaluate()

        moves = self.game.moves()
        maxval = -999
        for move in moves:
            self.game.make_move(move)
            val = self.alpha_beta(depth - 1)
            if val > maxval:
                maxval = val
            self.game.unmake_move()

        return maxval


##################
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
        # TODO: TESTING THIS. So right now, the rollout keeps generating same molecule
        #  after every reaction. I am assuming this is because most reactions we have now dont
        #  do anything so same product is formed after every reaction tghat is done. Needa
        #  test this with other reactions.
        current_rollout_state = self.state

        while not current_rollout_state.is_game_over():
            print("111", current_rollout_state.start)
            print(current_rollout_state.ctr)
            possible_moves = current_rollout_state.get_legal_actions()

            action = self.rollout_policy(possible_moves)
            temp = AlkeneReactions.partial_alkyne_hydrogenation
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
            print("NEW NODE CHOSEN FOR ROLLOUT", v.state.start) # The node chosen.
            if v.state.start == "CC=CCCCCC":
                print("Is da alkene")
                time.sleep(2)
            reward = v.rollout()
            v.backpropagate(reward)

        return self.best_child(c_param=0.)

class State:
    def __init__(self, state, end, ctr = 0):
        self.start = state
        self.end = end
        self.ctr = ctr # Denotes the no of times a reacn has been carried out

    def get_legal_actions(self):
        '''
        Modify according to your game or
        needs. Constructs a list of all
        possible actions from current state.
        Returns a list.
        '''
        # TODO: Maybe a better function here? Right now it returns all the rean.
        return AlkeneReactions.getAllReactions()

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
        # TODO: WE NEED A GOOD IS_GAME_OVER cuz the rollout never ends cuz eval <= 0.3 if unlucky
        evalu = self.evaluate(self.start, self.end)
        if self.ctr >= 30:
            return True

        if evalu == 1:
            return True

        if evalu <= 0.2:
            return False

    def game_result(self):
        '''
        Returns 1 or 0 or -1 depending
        on your state corresponding to win,
        tie or a loss.
        '''
        evalu = self.evaluate(self.start, self.end)
        if evalu == 1:
            return 1
        else:
            return -1



    def move(self, action):
        '''
        Changes the state of your
        board with a new value. Returns
        the new state after making a move.
        '''
        # self.start = action(self.start)[0]
        # self.ctr += 1
        # return self # TODO: Should this be a deepcopy or self. Not exactly sure

        newstart = action(self.start)[0] # This is making a deep copy and seems to work better
        newctr = self.ctr + 1
        return State(newstart, self.end, newctr)


def main(start, end):
    # TODO: Due to some reason, this seems to run reaaaallly slow. I an not sure why. Thinking
    #  cuz of all the new figures being created cuz bruh.
    # state = State("CC#CCCCCC", 'CCCCCC=O')
    state = State(start, end)
    root = MonteCarloTreeSearchNode(state=state)
    selected_node = root.best_action()
    x = selected_node
    print("FINAL", x.state.start)

    return x

def evaluate(smile1, smile2):
    ms = [Chem.MolFromSmiles(smile1), Chem.MolFromSmiles(smile2)]

    # for example
    fps = [Chem.RDKFingerprint(x) for x in ms]
    return DataStructs.FingerprintSimilarity(fps[0], fps[1])
