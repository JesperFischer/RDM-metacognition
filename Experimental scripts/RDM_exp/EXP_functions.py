"""
Created on Wed Feb 18 09:47:57 2026
@author: Siebe Everaerts

Code with functions for presenting random dot motion stimuli in 2 directions, recording responses
and giving feedback on the decision. Saves reaction time, response, accuracy, and
participant characteristics.

Use by pressing the "w" (or "z") key when dots are moving up and the "s" key when dots are moving down

"""

from psychopy import visual as vis
from psychopy import event
import numpy as np


################################################################################################################################ 
# Low Level functions
################################################################################################################################ 

## Function to check if the mouse is hovering over the slider bar area
def move_slider(left, right, up, SR, step_size, slider_pos):
    if left:
        slider_pos = max(0, slider_pos - step_size)  
    if right:
        slider_pos = min(1, slider_pos + step_size)
    if up:
        SR = slider_pos

    return slider_pos, SR


# Text for breaks in main experiment
def break_text_function(num_correct, tot_trials, mean_rt, per_correct, des_per_cor, des_mean_rt, blockN, n_blocks):
    points = 'You got ' + str(num_correct) + ' out of ' + str(tot_trials) + ' points this block!'
    speed = 'Your average reaction time was ' + str(f"{mean_rt:.2f}") + ' seconds'

    if per_correct < des_per_cor and mean_rt > des_mean_rt:
        feedback = 'Try to be faster and more accurate in the next block!'
    elif per_correct < des_per_cor and mean_rt <= des_mean_rt:
        feedback = 'Try to be more accurate in the next block!'
    elif per_correct >= des_per_cor and mean_rt > des_mean_rt:
        feedback = 'Try to be faster in the next block!'
    elif per_correct >= des_per_cor and mean_rt <= des_mean_rt:
        feedback = 'Good job! Try to be even faster and more accurate!'

    points_text = vis.TextStim(win, text = points, pos=(0, 100))
    speed_text = vis.TextStim(win, text = speed, pos = (0,50))
    break_text = vis.TextStim(win, text = "Take a short break before we continue with the next block (block " + str(blockN+1) + "/" + str(n_blocks) + ")", pos = (0, -100))
    space = vis.TextStim(win, text='Press space to continue', pos=(0, -150), height=20)
    feedback_text = vis.TextStim(win, text = feedback)
    return points_text, speed_text, break_text, space, feedback_text

# Z-score
def standard(val, mean,sd):
    Z = (val-mean)/sd
    return Z

# Psychometric weibull
def weibull_cdf(x, alpha, beta, gamma =0.5, lapse =0):
    return 1 - lapse - (1 - gamma - lapse) * (np.exp(-(x / alpha)**beta))

# Inverse weibull
def inv_weibull(y,alpha,beta):
    return alpha * (-np.log(2 * (1 - y)))**(1/beta)

def generate_valid_sequence(sequence, max_cons):
    
    def get_value(x):
        try:
            return x[0]      # If element is indexable (tuple/list/row), use first item
        except (TypeError, IndexError):
            return x         # If element is scalar, just use it
    
    while True:
        np.random.shuffle(sequence)   # Randomly reorder sequence
        
        counts = 1                    # Initialize consecutive count
        valid = True                  # Flag to track if sequence meets max repeat constraint
        
        for i in range(1, len(sequence)):
            if get_value(sequence[i]) == get_value(sequence[i-1]):  # Compare with previous element
                counts += 1            # Increment consecutive count
                if counts > max_cons:  # Check if repeats exceed limit
                    valid = False      # Sequence invalid, will reshuffle
                    break
            else:
                counts = 1             # Reset count when value changes
        
        if valid:
            return sequence.copy()     # Return a valid sequence without modifying original
        
