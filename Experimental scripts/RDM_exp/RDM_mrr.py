# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 09:47:57 2025
@author: Siebe Everaerts

Code for presenting random dot motion stimuli in 2 directions, recording responses
and giving feedback on the decision. Saves reaction time, response, accuracy, and
participant characteristics.

Use by pressing the "q" (or "a") key when dots are moving left and the "e" key when dots are moving right

"""
################################################################################################################################ 
# Packages
################################################################################################################################ 

# Set up experiment -----------------------------------------------------------
# Import modules
import random
import numpy as np
from psychopy import visual as vis
from psychopy import event, core, data
from psychopy.hardware import keyboard
from psychopy.data import QuestPlusHandler
import os
import instructions_1 as ins
import EXP_functions as exp
from scipy.stats import norm
import pylink
import sys
import serial
from pylink.eyelink import EyeLink
from EyeLinkCoreGraphicsPsychoPy import EyeLinkCoreGraphicsPsychoPy
import instructions_1 as ins


################################################################################################################################ 
# Set directory
################################################################################################################################ 
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)


################################################################################################################################ 
# Config
################################################################################################################################
# 
# ---------------------------------------------------------------
# Set up experiment parameters
# ---------------------------------------------------------------

keyboard_type = "qwerty"  # Keyboard layout ("qwerty" or "azerty")
pilot = 0                 # Pilot mode ON (1) or OFF (0)

# Staircase settings
SC_dotlife = 0            # Whether to run dotlife staircase (1 = yes, 0 = no)
SC_coherence = 0          # Whether to run coherence staircase (1 = yes, 0 = no)

# Training / instruction settings
tutorial = 0             # Run tutorial (1 = yes, 0 = no)
training_2 = 0            # Run second training block (1 = yes, 0 = no)
instructions = 0          # Show instructions (1 = yes, 0 = no)

fullscreen = 1            # Run experiment in fullscreen (1 = yes, 0 = no)
   

# ---------------------------------------------------------------
# Desired performance targets
# ---------------------------------------------------------------

des_per_cor = 0.6         # Desired percentage correct
des_mean_rt = 2           # Desired mean reaction time (seconds)
max_dur_conf = 3          # Maximum duration for confidence scale
step_size = 0.02          # Step size for confidence slider

# Accuracy values we want to simulate in the experiment
accVals = np.asarray([0.51, 0.6, 0.7, 0.8, 0.99])

# Number of repetitions for each coherence value in a sequence
repeats = [1, 2, 3, 2, 1]

# Maximum number of consecutive repetitions for up/down trials
max_cons = 4

# ---------------------------------------------------------------
# Trial characteristics
# ---------------------------------------------------------------

n_tutorial = 6 # Number of trials in second tutorial block

n_training_2 = 4 if pilot else 10  # Number of trials in second training block

n_seq = 2                           # Number of times to present the sequence
len_seq = 2 * sum(repeats)          # Length of one sequence (sum of repeats * 2)

# Total number of trials per testing block
n_trials = sum(repeats) if pilot else n_seq * len_seq

# Total number of blocks
n_blocks = 4 if pilot else 8

# Determine in which block other scale is presented
scale_Nblock = 1 if pilot else 4

# ---------------------------------------------------------------
# Hyperparameters for staircase dotlife
# ---------------------------------------------------------------

n_SC1 = 40           # Number of trials in dotlife staircase
coherence = 0.3      # Fixed coherence to find dotlife threshold performance

# ---------------------------------------------------------------
# Hyperparameters for staircase coherence
# ---------------------------------------------------------------

n_SC2 = 80           # Number of trials in coherence staircase


################################################################################################################################ 
# Setting up a data file
################################################################################################################################ 
# ---------------------------------------------------------------
# Participant information
# ---------------------------------------------------------------

if pilot:
    # Use fixed participant info for pilot mode
    sub = 0
    age = 30
    gender = 'Man'
    handedness = 'Right'
else:
    # Ask user to input participant info for real experiment
    sub = int(input("Subject number: "))
    age = int(input("Age: "))
    gender = input("Gender (Woman/Man/X): ")
    handedness = input("Handedness (Left/Right): ")

# Store participant info in a dictionary for easy access
info = {"sub": sub, "age": age, "gender": gender, "handedness": handedness}

# ---------------------------------------------------------------
# Setup file names for saving data
# ---------------------------------------------------------------

file_name = "Data/RDM_reportz_sub%d" % sub  # Main data file path
thisExp = data.ExperimentHandler(dataFileName=file_name,   extraInfo=info)


# ---------------------------------------------------------------
# Setup file names and path for physio
# ---------------------------------------------------------------
sys.path.append("Eyelink")
edf_remote_name = f"sqq_{sub}.EDF"
edf_local_path = os.path.join("Data", f"RDM_reportz_eyetrack_{sub}.EDF")


################################################################################################################################ 
# Psychopy objects
################################################################################################################################ 
# ---------------------------------------------------------------
# Set up Psychopy window
# ---------------------------------------------------------------
size = [1536, 960] if fullscreen else [600, 400]  # Window size depending on fullscreen mode
win = vis.Window(size=size,units='pix',color='grey',allowGUI=False,fullscr=fullscreen)

# Store window dimensions for later use
width = win.size[0]
height = win.size[1]

# ---------------------------------------------------------------
# Create core Psychopy objects
# ---------------------------------------------------------------
kb = keyboard.Keyboard()  # Keyboard object for responses
clock = core.Clock()      # Global clock for timing events

# ---------------------------------------------------------------
# Define response keys based on keyboard layout
# ---------------------------------------------------------------
if keyboard_type == "qwerty":
    choice_keys = ['w', 's', 'escape']  # up, down, escape
elif keyboard_type == "azerty":
    choice_keys = ['z', 's', 'escape']  # up, down, escape
else:
    raise TypeError('Unknown keyboard name')  # Error if layout unknown

# ---------------------------------------------------------------
# Create DotMotion stimulus
# Reference: https://www.psychopy.org/api/visual/dotstim.html
# ---------------------------------------------------------------
DotMotion = vis.DotStim(win, units='pix',nDots=120, fieldSize=300,fieldShape='circle',dotSize=7,
                        dotLife=15,speed=1,color='white',signalDots='same',    noiseDots='walk')

# ---------------------------------------------------------------
# Create slider for confidence or clarity ratings
# ---------------------------------------------------------------
slider = vis.Slider(win, name='slider',size=(400, 20), pos=(0, 0), units='pix',ticks=(0, 1), granularity=0.01,
    style=['rating'], color='white',font='HelveticaBold', flip=False)

slider.marker.color = "white"  # Slider marker color
slider.marker.size = 20        # Slider marker size

# Labels for confidence slider
slider_label_wrong = vis.TextStim(win, text="certainly wrong", pos=(-200, 30))
slider_label_right = vis.TextStim(win, text="certainly right", pos=(200, 30))

# Labels for clarity slider
slider_label_Nclear = vis.TextStim(win, text="not clear at all", pos=(-200, 30))
slider_label_clear = vis.TextStim(win, text="very clear", pos=(200, 30))

# Slider instructions
slider_instructions = vis.TextStim(win, text="How confident were you in your decision?", pos=(0, 75))
slider_instructions_dir = vis.TextStim(win, text="How clear was the DIRECTION of the dots?", pos=(0, 75))

# ---------------------------------------------------------------
# Text for no-response feedback
# ---------------------------------------------------------------
no_response_text = vis.TextStim(win,text="No response, try to be faster next trial!",color='white',height=30)

# ---------------------------------------------------------------
# Create fixation cross
# ---------------------------------------------------------------
fixation = vis.ShapeStim(win=win,vertices=((0, -15), (0, 15), (0, 0), (-15, 0), (15, 0)),  # Cross arms
    lineWidth=5,closeShape=False,lineColor='white')

# ---------------------------------------------------------------
# Space bar prompt for training / instructions
# ---------------------------------------------------------------
space = vis.TextStim(win, text='Press SPACE to continue', pos=(0, -300), height=30)


################################################################################################################################ 
# Set up Physio
################################################################################################################################ 

# ---------------------------------------------------------------
# Biopac
# ---------------------------------------------------------------
ser = serial.Serial('COM7', 115200, timeout=1)


# ---------------------------------------------------------------
# Eyelink
# ---------------------------------------------------------------
tracker = EyeLink("100.1.1.1")
tracker.openDataFile(edf_remote_name)

# EyeLink calibration graphics 
genv = EyeLinkCoreGraphicsPsychoPy(tracker, win)
genv.setCalibrationColors((-1, -1, -1), win.color)
genv.setTargetType('picture')
genv.setPictureTarget(os.path.join(dname, 'Eyelink_cal', 'fixTarget.bmp'))
pylink.openGraphicsEx(genv)

# EyeLink screen and calibration settings 
tracker.sendCommand(f"screen_pixel_coords = 0 0 {width - 1} {height - 1}")
tracker.sendMessage(f"DISPLAY_COORDS = 0 0 {width - 1} {height - 1}")
tracker.sendCommand('enable_automatic_calibration=YES')
tracker.sendCommand('automatic_calibration_pacing=500')
           
# ---------------------------------------------------------------
# Eyelink calibration
# ---------------------------------------------------------------
# 
tracker.doTrackerSetup()

################################################################################################################################ 
# tutorial block
################################################################################################################################ 
ins.Intro(win); core.wait(1); event.waitKeys(keyList=['space'])

block = 0
TrialType = "tutorial"
coherence = 0.9

acc = [0] * n_tutorial

if tutorial: 
    # Balanced left/right sequence for training
    sequence = n_tutorial // 2 * [0] + n_tutorial // 2 * [1]
    condition_direction = exp.generate_valid_sequence(sequence, max_cons)  # Avoid too many repeats

    # Loop through all trials
    for trial in range(n_tutorial):
        
        print("trial start")

        # Map 0/1 to motion direction and angle
        mapping = {0: ('up', 90), 1: ('down', 270)}
        correct, direction = mapping[condition_direction[trial]]    

        # Prepare for trial
        resp = None          # Empty variable for response
        event.clearEvents()  # Clear any lingering key presses

        # Assign stimulus parameters
        DotMotion.coherence = coherence   # constant coherence
        DotMotion.dir = direction         # Set motion direction

        # ---------------------------------------------------------------
        # Present stimulus and wait for response
        # ---------------------------------------------------------------
        while not resp:
            fixation.draw()     # Draw fixation cross
            DotMotion.draw()    # Draw random dot stimulus
            win.flip()
            resp = event.getKeys(keyList=choice_keys)  # Collect response

        # ---------------------------------------------------------------
        # Calculate accuracy
        # ---------------------------------------------------------------
        correct_key = choice_keys[0] if correct == "up" else choice_keys[1]
        if resp:
            is_correct = (resp[0] == correct_key)
            ACC = int(is_correct)
            acc[trial] = ACC

            # Feedback text
            if is_correct:
                print("Decision was correct")
                FB_text = "Correct!"
            else:
                print("Decision was incorrect")
                FB_text = "Wrong"
        else:
            ACC = 0  # No response = incorrect

        # ---------------------------------------------------------------
        # Allow escape to exit experiment
        # ---------------------------------------------------------------
        if resp == ['escape']:
            print('Participant pressed escape')
            thisExp.saveAsWideText(file_name + '.csv', delim=',')
            win.close()
            core.quit()

        # ---------------------------------------------------------------
        # Display feedback
        # ---------------------------------------------------------------
        feedback = vis.TextStim(win, text=FB_text, color="white", height=40)
        feedback.draw()
        win.flip()
        core.wait(0.5)

        # ---------------------------------------------------------------
        # Convert key press to "up"/"down" label
        # ---------------------------------------------------------------
        key_to_label = {choice_keys[0]: "up", choice_keys[1]: "down"}
        if resp:
            resp = key_to_label[resp[0]]

        # ---------------------------------------------------------------
        # Save trial data to experiment handler
        # ---------------------------------------------------------------
        thisExp.addData("block", block)
        thisExp.addData("Trialtype", TrialType)
        thisExp.addData("withinblocktrial", trial)
        thisExp.addData("resp", resp)
        thisExp.addData("cor", ACC)
        thisExp.addData("dots direction", direction)
        thisExp.addData("cor_resp", correct)
        thisExp.addData("coherence", DotMotion.coherence)
        thisExp.nextEntry()


ins.Tutorial(win); core.wait(1); event.waitKeys(keyList=['space'])


################################################################################################################################ 
# Staircase Dotlife
################################################################################################################################
# ---------------------------------------------------------------
# Hyperparameters & Staircase setup
# ---------------------------------------------------------------

TrialType = "SC dotlife"  # Label for this trial type
block =+ 1                # Current block number

# Create a balanced sequence for left/right directions
sequence = n_SC1 // 2 * [0] + n_SC1 // 2 * [1]

# Priors for QUEST+ staircase (normalize PDFs)
threshold_prior = norm.pdf(np.linspace(1, 20, 50), loc=10, scale=5)
threshold_prior = threshold_prior / threshold_prior.sum()

slope_prior = norm.pdf(np.linspace(1, 10, 10), loc=5, scale=2)
slope_prior = slope_prior / slope_prior.sum()

# Initialize QUEST+ staircase
SC = QuestPlusHandler(nTrials=n_SC1,intensityVals=list(range(1, 20, 2)),thresholdVals=np.linspace(1, 20, 50),slopeVals=np.linspace(1, 10, 10),lowerAsymptoteVals=0.5,
    lapseRateVals=0.03,responseVals=[1, 0],prior={"threshold": threshold_prior, "slope": slope_prior},psychometricFunc="weibull",startIntensity=11,
    stimScale="linear",stimSelectionMethod="minEntropy",paramEstimationMethod="mean")

# ---------------------------------------------------------------
# Run staircase trials if SC_dotlife is True
# ---------------------------------------------------------------
if SC_dotlife:
    
    # Balanced left/right sequence for training
    sequence = n_SC1 // 2 * [0] + n_SC1 // 2 * [1]
    condition_direction = exp.generate_valid_sequence(sequence, max_cons)  # Avoid too many repeats

    # Initialize empty lists for accuracy and reaction times
    acc = [0] * n_SC1
    rt = [0] * n_SC1

    # Loop through all trials
    for trial in range(n_SC1):
        
        print("trial start")
        event.clearEvents()  # Clear keyboard events at start
        
        # Get staircase value for this trial
        stim = SC.next()
        print("extracting SC value successful")
        
        # Map 0/1 to motion direction and angle
        mapping = {0: ('up', 90), 1: ('down', 270)}
        correct, direction = mapping[condition_direction[trial]]    

        # Prepare for trial
        resp = None          # Empty variable for response
        event.clearEvents()  # Clear any lingering key presses

        # Assign stimulus parameters
        DotMotion.coherence = coherence   # constant coherence
        DotMotion.dotLife = stim          # dotLife set by staircase
        print(DotMotion.dotLife)
        DotMotion.dir = direction         # Set motion direction

        # ---------------------------------------------------------------
        # Present stimulus and wait for response
        # ---------------------------------------------------------------
        T_stimulus_start = clock.getTime()
        while not resp:
            fixation.draw()     # Draw fixation cross
            DotMotion.draw()    # Draw random dot stimulus
            win.flip()
            resp = event.getKeys(keyList=choice_keys)  # Collect response

            # Timeout: skip trial if no response
            if clock.getTime() - T_stimulus_start >= des_mean_rt:
                print("No response within 2 s, skipping trial")
                FB_text = "No response"
                FB_col = "white"
                break
        
        # ---------------------------------------------------------------
        # Record reaction time
        # ---------------------------------------------------------------
        if resp:
            T_stimulus_stop = clock.getTime()
            RTdec = T_stimulus_stop - T_stimulus_start
            rt[trial] = RTdec
            print("response given")
            print("Reaction time is:", RTdec)
        else:
            RTdec = np.nan
            rt[trial] = RTdec

        # ---------------------------------------------------------------
        # Calculate accuracy
        # ---------------------------------------------------------------
        correct_key = choice_keys[0] if correct == "up" else choice_keys[1]
        if resp:
            is_correct = (resp[0] == correct_key)
            ACC = int(is_correct)
            acc[trial] = ACC

            # Feedback text
            if is_correct:
                print("Decision was correct")
                FB_text = "Correct!"
            else:
                print("Decision was incorrect")
                FB_text = "Wrong"
        else:
            ACC = 0  # No response = incorrect

        # ---------------------------------------------------------------
        # Allow escape to exit experiment
        # ---------------------------------------------------------------
        if resp == ['escape']:
            print('Participant pressed escape')
            thisExp.saveAsWideText(file_name + '.csv', delim=',')
            win.close()
            core.quit()

        # ---------------------------------------------------------------
        # Update staircase
        # ---------------------------------------------------------------
        SC.addResponse(ACC)
        print("update staircase = OK")

        # ---------------------------------------------------------------
        # Display feedback
        # ---------------------------------------------------------------
        feedback = vis.TextStim(win, text=FB_text, color="white", height=40)
        feedback.draw()
        win.flip()
        core.wait(0.5 if resp else 1)  # Show feedback slightly longer if no response

        # ---------------------------------------------------------------
        # Convert key press to "up"/"down" label
        # ---------------------------------------------------------------
        key_to_label = {choice_keys[0]: "up", choice_keys[1]: "down"}
        if resp:
            resp = key_to_label[resp[0]]

        # ---------------------------------------------------------------
        # Save trial data to experiment handler
        # ---------------------------------------------------------------
        thisExp.addData("block", block)
        thisExp.addData("Trialtype", TrialType)
        thisExp.addData("withinblocktrial", trial)
        thisExp.addData("RTdec", RTdec)
        thisExp.addData("resp", resp)
        thisExp.addData("cor", ACC)
        thisExp.addData("dots direction", direction)
        thisExp.addData("cor_resp", correct)
        thisExp.addData("coherence", DotMotion.coherence)
        thisExp.addData("dotlife", DotMotion.dotLife)
        thisExp.nextEntry()
     

################################################################################################################################ 
# Staircase Coherence
################################################################################################################################ 
# ---------------------------------------------------------------
# Hyperparameters & Staircase setup for Coherence
# ---------------------------------------------------------------

TrialType = "SC coherence"  # Label for this trial type
block += 1
n_SC2 = 80  # Number of trials in this staircase

# Balanced left/right sequence
sequence = n_SC2 // 2 * [0] + n_SC2 // 2 * [1]

# Set dotLife to the previously estimated threshold (task difficulty fixed)
dotLife = SC.paramEstimate["threshold"]
print("Dotlife threshold: ", dotLife)

# Priors for QUEST+ staircase (normalize PDFs)
threshold_prior = norm.pdf(np.arange(0, 1, 0.02), loc=0.3, scale=0.2)
threshold_prior = threshold_prior / threshold_prior.sum()

slope_prior = norm.pdf(np.arange(0.5, 10.1, 0.5), loc=4, scale=2)
slope_prior = slope_prior / slope_prior.sum()

# Initialize QUEST+ staircase for coherence
SC = QuestPlusHandler(nTrials=n_SC2,intensityVals=np.arange(0, 1, 0.02),thresholdVals=np.arange(0, 1, 0.02),slopeVals=np.arange(0.5, 10.1, 0.5),lowerAsymptoteVals=0.5,lapseRateVals=0,
    responseVals=[1, 0],prior={"threshold": threshold_prior, "slope": slope_prior},psychometricFunc="weibull",startIntensity=0.3,stimScale="linear",stimSelectionMethod="minEntropy",paramEstimationMethod="mean")



# ---------------------------------------------------------------
# Run SC coherence staircase trials
# ---------------------------------------------------------------
if SC_coherence:

    # Create balanced left/right directions
    condition_direction = exp.generate_valid_sequence(sequence, max_cons)
    
    # Empty lists for accuracy and reaction times
    acc = [0] * n_SC2
    rt = [0] * n_SC2

    # Loop through all trials
    for trial in range(n_SC2):

        # Get staircase value for this trial (coherence)
        stim = SC.next()
        print("Extracting SC value successful")

        # Map 0/1 to motion direction and angle
        mapping = {0: ('up', 90), 1: ('down', 270)}
        correct, direction = mapping[condition_direction[trial]]

        # Prepare trial
        resp = None
        event.clearEvents()  # Clear previous keyboard events

        # Assign stimulus parameters
        DotMotion.coherence = stim   # Coherence set by staircase
        DotMotion.dotLife = dotLife  # DotLife fixed
        DotMotion.dir = direction
        print(stim)    # Print fixed dotLife for verification
        print("Presenting stimulus...")  # Print trial start

        # ---------------------------------------------------------------
        # Present stimulus and wait for response
        # ---------------------------------------------------------------
        T_stimulus_start = clock.getTime()
        while not resp:
            fixation.draw()
            DotMotion.draw()
            win.flip()
            resp = event.getKeys(keyList=choice_keys)

            # Timeout: skip trial if no response
            if clock.getTime() - T_stimulus_start >= des_mean_rt:
                print("No response within 2 s, skipping trial")
                FB_text = "No response"
                FB_col = "white"
                break

        # ---------------------------------------------------------------
        # Record reaction time
        # ---------------------------------------------------------------
        if resp:
            T_stimulus_stop = clock.getTime()
            RTdec = T_stimulus_stop - T_stimulus_start
            rt[trial] = RTdec
            print("Response given")
            print("Reaction time is:", RTdec)
        else:
            RTdec = np.nan
            rt[trial] = RTdec

        # ---------------------------------------------------------------
        # Calculate accuracy
        # ---------------------------------------------------------------
        correct_key = choice_keys[0] if correct == "up" else choice_keys[1]
        if resp:
            is_correct = (resp[0] == correct_key)
            ACC = int(is_correct)
            acc[trial] = ACC

            # Feedback text
            if is_correct:
                print("Decision was correct")
                FB_text = "Correct!"
            else:
                print("Decision was incorrect")
                FB_text = "Wrong"
        else:
            ACC = 0  # No response = incorrect

        # ---------------------------------------------------------------
        # Allow escape to exit experiment
        # ---------------------------------------------------------------
        if resp == ['escape']:
            print('Participant pressed escape')
            thisExp.saveAsWideText(file_name + '.csv', delim=',')
            win.close()
            core.quit()

        # ---------------------------------------------------------------
        # Update staircase
        # ---------------------------------------------------------------
        SC.addResponse(ACC)
        print("Update staircase = OK")

        # ---------------------------------------------------------------
        # Display feedback
        # ---------------------------------------------------------------
        feedback = vis.TextStim(win, text=FB_text, color="white", height=40)
        feedback.draw()
        win.flip()
        core.wait(0.5 if resp else 1)

        # ---------------------------------------------------------------
        # Convert key press to "up"/"down" label
        # ---------------------------------------------------------------
        key_to_label = {choice_keys[0]: "up", choice_keys[1]: "down"}
        if resp:
            resp = key_to_label[resp[0]]

        # ---------------------------------------------------------------
        # Save trial data
        # ---------------------------------------------------------------
        thisExp.addData("block", block)
        thisExp.addData("Trialtype", TrialType)
        thisExp.addData("withinblocktrial", trial)
        thisExp.addData("RTdec", RTdec)
        thisExp.addData("resp", resp)
        thisExp.addData("cor", ACC)
        thisExp.addData("dots direction", direction)
        thisExp.addData("cor_resp", correct)
        thisExp.addData("coherence", DotMotion.coherence)
        thisExp.addData("dotlife", DotMotion.dotLife)
        thisExp.nextEntry()

# ---------------------------------------------------------------
# Extract final slope and threshold
# ---------------------------------------------------------------
print(SC.paramEstimate)
threshold_coh = SC.paramEstimate["threshold"]
slope_coh = SC.paramEstimate["slope"]

# Calculate coherence values for the real experiment
coherenceVals = np.minimum(exp.inv_weibull(accVals, threshold_coh, slope_coh), 1)
print(coherenceVals)

################################################################################################################################ 
# Instructions
################################################################################################################################ 

# ---------------------------------------------------------------
# Show instruction screens
# ---------------------------------------------------------------
if instructions:

    # Sequential instruction pages (advance with space)
    ins.Main1(win); core.wait(1); event.waitKeys(keyList=['space'])
    ins.Main2(win); core.wait(1); event.waitKeys(keyList=['space'])
    ins.Main3(win); core.wait(1); event.waitKeys(keyList=['space'])
    ins.Main4(win); core.wait(1); event.waitKeys(keyList=['space'])


    # ---------------------------------------------------------------
    # Create slider practice instructions
    # ---------------------------------------------------------------

    # Top instruction: explain how to move the slider
    ins_scale = vis.TextStim(win=win,text="Try playing around with the scale. Press ← and → to move the cursor around.",
        pos=(0, 300),color='white',height=38,wrapWidth=1400,alignText='center')

    # Bottom instruction: explain how to confirm
    ins_scale_1 = vis.TextStim(win=win,text="Press ↑ whenever you are ready to continue.",
        pos=(0, -300),color='white',height=38,wrapWidth=1400,alignText='center')
    
    # ---------------------------------------------------------------
    # Prepare slider practice screen
    # ---------------------------------------------------------------

    event.clearEvents()      # Clear old keyboard events
    kb.clearEvents()         # Clear keyboard state

    slider.reset()           # Reset slider to initial state
    slider.markerPos = 0.5   # Start marker in middle of scale

    # Draw instructions and slider
    ins_scale.draw()
    ins_scale_1.draw()
    slider.draw()
    win.flip()


    # ---------------------------------------------------------------
    # Allow participant to practice moving slider
    # Loop until confirmation key (up arrow) is pressed
    # ---------------------------------------------------------------

    SR = None                # Will store confirmation signal
    held_keys = []           # (Optional) track held keys if used in move function

    while SR is None:

        # -----------------------------------------------------------
        # Read keyboard state (continuous polling)
        # -----------------------------------------------------------
        left, right, up, escape = kb.getState(['left', 'right', 'up', 'escape'])


        # -----------------------------------------------------------
        # Update slider position based on key input
        # -----------------------------------------------------------
        slider_pos = slider.markerPos

        slider_pract_value, SR = exp.move_slider(left,right,up,SR,step_size,slider_pos)

        # Apply updated slider position
        slider.markerPos = slider_pract_value

        # -----------------------------------------------------------
        # Redraw screen each frame
        # -----------------------------------------------------------
        slider.draw();ins_scale.draw();ins_scale_1.draw();win.flip()

    # ---------------------------------------------------------------
    # Continue to final instruction page
    # ---------------------------------------------------------------
    ins.Main5(win);core.wait(1);event.waitKeys(keyList=['space'])
################################################################################################################################ 
# Training 2
################################################################################################################################
# ---------------------------------------------------------------
# Start training block with confidence scale
# ---------------------------------------------------------------

block += 1
TrialType = "Training - scale"


if training_2:

    # -----------------------------------------------------------
    # Generate stimulus parameters for training trials
    # -----------------------------------------------------------

    # Draw coherence values and randomize order
    coherence = np.linspace(0.1, 1.0, n_training_2)
    random.shuffle(coherence)

    # Create balanced up/down sequence
    sequence = n_training_2 // 2 * [0] + n_training_2 // 2 * [1]

    # Generate direction order with max repetition constraint
    condition_direction = exp.generate_valid_sequence(sequence, max_cons)

    # Prepare storage for performance measures
    acc = [0] * n_training_2
    rt = [0] * n_training_2


    # -----------------------------------------------------------
    # Loop through all training trials
    # -----------------------------------------------------------
    for trial in range(n_training_2):

        # -------------------------------------------------------
        # Determine stimulus direction for this trial
        # -------------------------------------------------------
        mapping = {0: ('up', 90), 1: ('down', 270)}
        correct, direction = mapping[condition_direction[trial]]  


        # -------------------------------------------------------
        # Prepare stimulus presentation
        # -------------------------------------------------------
        resp = None
        event.clearEvents()

        DotMotion.coherence = coherence[trial]
        DotMotion.dotLife = dotLife
        DotMotion.dir = direction

        print("Presenting stimulus...")


        # -------------------------------------------------------
        # Present stimulus and wait for response (decision phase)
        # -------------------------------------------------------
        T_stimulus_start = clock.getTime()

        while not resp:
            fixation.draw()
            DotMotion.draw()
            win.flip()

            resp = event.getKeys(keyList=choice_keys)

            # Timeout → show feedback and skip trial
            if clock.getTime() - T_stimulus_start >= des_mean_rt:
                print("No response within 2 s, skipping trial")
                miss_text = vis.TextStim(
                    win,
                    text="No response, try to be faster next trial!",
                    height=30
                )
                miss_text.draw()
                win.flip()
                core.wait(1)
                break


        # -------------------------------------------------------
        # Record reaction time
        # -------------------------------------------------------
        if resp:
            T_stimulus_stop = clock.getTime()
            RTdec = T_stimulus_stop - T_stimulus_start
            rt[trial] = RTdec
            print("participant responded")
            print("Reaction time is:", RTdec)
        else:
            RTdec = np.nan
            rt[trial] = RTdec


        # -------------------------------------------------------
        # Compute decision accuracy
        # -------------------------------------------------------
        correct_key = choice_keys[0] if correct == "up" else choice_keys[1]

        if resp:
            is_correct = (resp[0] == correct_key)
            ACC = int(is_correct)
            acc[trial] = ACC
        else:
            ACC = 0


        # -------------------------------------------------------
        # Allow escape to abort experiment safely
        # -------------------------------------------------------
        if resp == ['escape']:
            print('Participant pressed escape')
            thisExp.saveAsWideText(file_name + '.csv', delim=',')
            win.close()
            core.quit()


        # -------------------------------------------------------
        # Inter-stimulus interval before confidence rating
        # -------------------------------------------------------
        fixation.draw()
        win.flip()
        core.wait(1)


        # -------------------------------------------------------
        # Prepare confidence slider
        # -------------------------------------------------------
        kb.clearEvents()
        slider.reset()

        # Randomize initial marker position (reduces anchoring bias)
        slider.markerPos = np.random.uniform(0.25, 0.75)

        print("confidence start")


        # -------------------------------------------------------
        # Confidence rating phase (only if response was given)
        # -------------------------------------------------------
        if resp:

            slider.draw()
            slider_instructions.draw()
            slider_label_wrong.draw()
            slider_label_right.draw()
            win.flip()

            SR = None
            held_keys = []

            # Loop until participant confirms rating
            while SR is None:

                # Read continuous key state
                left, right, up, escape = kb.getState(['left', 'right', 'up', 'escape'])

                # Emergency exit
                if escape:
                    print('Participant pressed escape')
                    thisExp.saveAsWideText(file_name + '.csv', delim=',')
                    win.close()
                    core.quit()

                # Update slider position
                slider_pos = slider.markerPos
                slider_pract_value, SR = exp.move_slider(
                    left,
                    right,
                    up,
                    SR,
                    step_size,
                    slider_pos
                )

                slider.markerPos = slider_pract_value

                # Redraw rating screen
                slider.draw()
                slider_instructions.draw()
                slider_label_wrong.draw()
                slider_label_right.draw()
                win.flip()

            print("participant answered")
            print("Reported confidence = ", SR)

        else:
            # No decision → no confidence rating
            RTrating = None
            SR = None


        # -------------------------------------------------------
        # Inter-trial interval
        # -------------------------------------------------------
        fixation.draw()
        win.flip()
        core.wait(1)


        # -------------------------------------------------------
        # Convert response key to label
        # -------------------------------------------------------
        key_to_label = {choice_keys[0]: "up", choice_keys[1]: "down"}

        if resp:
            resp = key_to_label[resp[0]]


        # -------------------------------------------------------
        # Save trial data
        # -------------------------------------------------------
        thisExp.addData("block", block)
        thisExp.addData("Trialtype", TrialType)
        thisExp.addData("withinblocktrial", trial)
        thisExp.addData("RTdec", RTdec)
        thisExp.addData("resp", resp)
        thisExp.addData("cor", ACC)
        thisExp.addData("dots direction", direction)
        thisExp.addData("cor_resp", correct)
        thisExp.addData("SR_conf", SR)
        thisExp.addData("coherence", coherence[trial])
        thisExp.addData("dotlife", dotLife)
        thisExp.nextEntry()


# ---------------------------------------------------------------
# End of training block → show next instruction screen
# ---------------------------------------------------------------
ins.Main6(win); core.wait(1); event.waitKeys(keyList=['space'])

################################################################################################################################ 
# Variables for main experiment
################################################################################################################################
TrialType = "Main"  # Define trial type for logging


# ---------------------------------------------------------------
# Prepare stimuli for first block
# ---------------------------------------------------------------

# Create repeated coherence values for each trial
sequence = np.repeat(coherenceVals, repeats)
sequence = np.concatenate([sequence, sequence])

# Create balanced up/down directions
direction = [1] * sum(repeats) + [0] * sum(repeats)

# Combine direction and coherence into stimulus info
stim_info_1 = list(zip(direction, sequence))
stim_info = []

# Generate valid sequences (avoid too many repeats) for multiple repetitions
for _ in range(n_seq):
    valid_seq = exp.generate_valid_sequence(stim_info_1, max_cons)
    stim_info.extend(valid_seq)


# ---------------------------------------------------------------
# Determine inter-trial waiting times
# ---------------------------------------------------------------

# Inter-trial interval drawn from a log-normal style distribution
inter_trial = np.exp(np.random.normal(loc=1, scale=0.1, size=n_trials))

# Interval for confidence rating based on previous studies
manipulation = np.random.uniform(low=des_mean_rt, high=5, size=n_trials)


# ---------------------------------------------------------------
# Initialize performance vectors and counters
# ---------------------------------------------------------------
acc = [0] * n_trials
rt = [0] * n_trials
trialN = 0
blockN = 0


################################################################################################################################ 
# Main experiment
################################################################################################################################ 

# ---------------------------------------------------------------
# Loop through all trials across all blocks
# ---------------------------------------------------------------
for eachTrial in range(n_trials * n_blocks):

    # -----------------------------------------------------------
    # Assign stimulus parameters for this trial
    # -----------------------------------------------------------
    mapping = {0: ('up', 90), 1: ('down', 270)}
    correct, direction = mapping[stim_info[trialN][0]]

    # Save stimulus start time
    T_stimulus_start = clock.getTime()

    resp = None
    event.clearEvents()

    DotMotion.coherence = stim_info[trialN][1]
    DotMotion.dotLife = dotLife
    DotMotion.dir = direction

    print("presenting stimulus ...")

    # start the eyetracker:
    tracker.startRecording(1, 1, 1, 1)
    tracker.sendMessage(f"start_trialID_{trialN}_Block_{blockN + 2}")

    # biopack
    ser.write(str.encode('01'))
    core.wait(0.1)
    ser.write(str.encode('00')) 
    
    tracker.sendMessage("start_stimulus")  # Optional: timestamp visual onset

    # -----------------------------------------------------------
    # Present stimulus until response or timeout
    # -----------------------------------------------------------
    while not resp:
        fixation.draw()
        DotMotion.draw()
        win.flip()

        resp = event.getKeys(keyList=choice_keys)

        # Timeout: show message and skip trial
        if clock.getTime() - T_stimulus_start >= des_mean_rt:
            print("No response within 1.5 s, skipping trial")
            miss_text = vis.TextStim(
                win,
                text="No response, try to be faster next trial!",
                height=30
            )
            miss_text.draw()
            win.flip()
            core.wait(1)
            break
    
    # send message to tracker
    tracker.sendMessage("end_stimulus")  

    # -----------------------------------------------------------
    # Record reaction time
    # -----------------------------------------------------------
    if resp:
        T_stimulus_stop = clock.getTime()
        RTdec = T_stimulus_stop - T_stimulus_start
        rt[trialN] = RTdec
        print("participant responded")
        print("Reaction time is:", RTdec)
    else:
        rt[trialN] = np.nan


    # -----------------------------------------------------------
    # Compute accuracy
    # -----------------------------------------------------------
    correct_key = choice_keys[0] if correct == "up" else choice_keys[1]

    if resp:
        is_correct = (resp[0] == correct_key)
        ACC = int(is_correct)
        acc[trialN] = ACC
    else:
        ACC = 0
            
    # -----------------------------------------------------------
    # Emergency exit: save and quit if escape is pressed
    # -----------------------------------------------------------
    if resp == ['escape']:
        print('Participant pressed escape')
        thisExp.saveAsWideText(file_name + '.csv', delim=',')

        # save eyelink EDF from tracker to local Data folder
        tracker.setOfflineMode()
        tracker.closeDataFile()
        try:
            print("Receiving EDF from EyeLink...")
            tracker.receiveDataFile(edf_remote_name, edf_local_path)
            print(f"EDF saved to {edf_local_path}")
        except RuntimeError as e:
            print("Error transferring EDF:", e)

        tracker.close()

        win.close()
        core.quit()

    # -----------------------------------------------------------
    # Confidence rating phase
    # -----------------------------------------------------------
    print("confidence start")
    interval = manipulation[trialN]

    fixation.draw()
    win.flip()

    if resp:
        # Wait for confidence interval minus decision RT
        core.wait(interval - rt[trialN] - 0.05)
        ser.write(str.encode('01'))
        core.wait(0.05)
        ser.write(str.encode('00')) # turn off all 8

        # Initialize confidence slider
        T_rating_start = clock.getTime()
        slider.draw()
        if blockN == scale_Nblock:
            slider_instructions_dir.draw();slider_label_Nclear.draw();slider_label_clear.draw()
        else:
            slider_instructions.draw();slider_label_wrong.draw();slider_label_right.draw()

        win.flip()

        SR = None
        slider.reset()
        start_conf = np.random.uniform(0.25,0.75)
        slider.markerPos = start_conf

        # Loop until participant confirms rating or timeout
        while SR is None:
            elapsed_time = clock.getTime() - T_rating_start
            if elapsed_time >= max_dur_conf:
                no_response_text.draw()
                win.flip()
                core.wait(1)
                SR = None
                RTrating = None
                break

            # Get key states
            left, right, up, escape = kb.getState(['left', 'right', 'up', 'escape'])

                        # Escape check
            if escape:
                print('Participant pressed escape')
                thisExp.saveAsWideText(file_name + '.csv', delim=',')

                tracker.setOfflineMode()
                tracker.closeDataFile()
                try:
                    print("Receiving EDF from EyeLink...")
                    tracker.receiveDataFile(edf_remote_name, edf_local_path)
                    print(f"EDF saved to {edf_local_path}")
                except RuntimeError as e:
                    print("Error transferring EDF:", e)
                tracker.close()
                win.close()  
                core.quit()


            # -----------------------------------------------------------
            # Emergency exit: save and quit if escape is pressed
            # -----------------------------------------------------------
            if escape:  
                print('Participant pressed escape')
                thisExp.saveAsWideText(file_name + '.csv', delim=',') 
                                ## save eyelink EDF from tracker to local Data folder
                #tracker.setOfflineMode()
                #tracker.closeDataFile()
                try:
                    print("Receiving EDF from EyeLink...")
                    #tracker.receiveDataFile(edf_remote_name, edf_local_path)
                    print(f"EDF saved to {edf_local_path}")
                except RuntimeError as e:
                    print("Error transferring EDF:", e)
                #tracker.close()
                win.close()  
                core.quit()

            # Update slider position
            slider_pos = slider.markerPos
            slider_pract_value, SR = exp.move_slider(left, right, up, SR, step_size, slider_pos)
            slider.markerPos = slider_pract_value

            # Redraw slider and instructions
            slider.draw()
            if blockN == scale_Nblock:
                slider_instructions_dir.draw()
                slider_label_Nclear.draw()
                slider_label_clear.draw()
                scale = 'control'
            else:
                slider_instructions.draw()
                slider_label_wrong.draw()
                slider_label_right.draw()
                scale = 'conf'
            win.flip()

            T_rating_stop = clock.getTime()
            RTrating = T_rating_stop - T_rating_start

    else:
        RTrating = None
        SR = None
        interval = None

    
    # -----------------------------------------------------------
    # Post-confidence ITI
    # -----------------------------------------------------------
    print("participant responded")
    core.wait(0.1)
    fixation.draw()
    win.flip()
    core.wait(inter_trial[trialN])
    
    core.wait(0.1)
    tracker.sendMessage(f"end_trialID_{trialN}_Block_{blockN + 2}")
    ser.write(str.encode('01'))
    core.wait(0.015)
    ser.write(str.encode('00'))
    tracker.stopRecording()
        

    # -----------------------------------------------------------
    # Convert key press to label for saving
    # -----------------------------------------------------------
    key_to_label = {choice_keys[0]: "up", choice_keys[1]: "down"}
    if resp:
        resp = key_to_label[resp[0]]
    else:
        scale = None


    # -----------------------------------------------------------
    # Save trial data
    # -----------------------------------------------------------
    thisExp.addData("block", blockN + 2)
    thisExp.addData("Trialtype", TrialType)
    thisExp.addData("withinblocktrial", trialN)
    thisExp.addData("RTdec", RTdec)
    thisExp.addData("resp", resp)
    thisExp.addData("cor", ACC)
    thisExp.addData("dots direction", direction)
    thisExp.addData("cor_resp", correct)
    thisExp.addData("interval", interval)
    thisExp.addData("interTrial.interval", inter_trial[trialN])
    thisExp.addData("scale", scale)
    thisExp.addData("start_conf", start_conf)
    thisExp.addData("SR_conf", SR)
    thisExp.addData("RTrating", RTrating)
    thisExp.addData("coherence", DotMotion.coherence)
    thisExp.addData("dotlife", dotLife)
    thisExp.nextEntry()
      # -----------------------------------------------------------
    # Update trial and block counters
    # -----------------------------------------------------------
    trialN += 1

    # Move to next block if finished current block
    if blockN < n_blocks - 1 and trialN == n_trials:
        blockN += 1
        trialN = 0

        # Recompute inter-trial and confidence intervals
        inter_trial = np.exp(np.random.normal(loc=1, scale=0.1, size=n_trials))
        manipulation = np.random.uniform(low=des_mean_rt, high=5, size=n_trials)

        # Prepare stimuli for next block
        sequence = np.repeat(coherenceVals, repeats)
        sequence = np.concatenate([sequence, sequence])
        direction = [1] * 9 + [0] * 9

        stim_info_1 = list(zip(direction, sequence))
        stim_info = []
        for _ in range(n_seq):
            valid_seq = exp.generate_valid_sequence(stim_info_1)
            stim_info.extend(valid_seq)

        # Performance summary for current block
        num_correct = sum(acc)
        tot_trials = len(acc)
        per_correct = sum(acc) / len(acc)
        mean_rt = np.nanmean(rt)

        # Reset performance vectors
        acc = [0] * n_trials
        rt = [0] * n_trials
                
        # -------------------------------------------------------
        # Show break and feedback
        # -------------------------------------------------------
        points_text, speed_text, break_text, space, feedback_text = exp.break_text_function(num_correct, tot_trials, mean_rt, per_correct,des_per_cor, des_mean_rt, blockN, n_blocks)
        points_text.draw(); speed_text.draw(); feedback_text.draw(); break_text.draw(); win.flip(); core.wait(5)
        points_text.draw(); speed_text.draw(); feedback_text.draw(); break_text.draw(); space.draw(); win.flip(); event.waitKeys(keyList=['space'])

        # Extra instructions depending on block
        if blockN == scale_Nblock:
            ins.ExtraScale(win); core.wait(1); event.waitKeys(keyList=['lctrl'])
        elif blockN == scale_Nblock + 1:
            ins.Main7(win); core.wait(1); event.waitKeys(keyList=['lctrl'])

        # Show recaliberation screen
        recalibration_text = vis.TextStim(win, text = "Awaiting recalibration... \n Call the experimenter.", pos = (0,0), height = 40, wrapWidth = 1200)
        recalibration_text.draw(); win.flip(); event.waitKeys(keyList=['lctrl']) 
        
        # recalibration
        tracker.sendMessage("recalibration_start")
        tracker.doTrackerSetup()  # recalibrate eyetracker at break
        tracker.sendMessage("recalibration_end")

        next_block = vis.TextStim(win, text = "Please put your hands back on the keyboard. \n Get ready to restart :)", color = "white", height=40)
        next_block.draw(); win.flip(); core.wait(5)
        
################################################################################################################################ 
# End of experiment
################################################################################################################################ 

# Thank participant and give instructions to close experiment
break_text = vis.TextStim(win,text="This is the end of the experiment. \n Thank you very much for your participation!")
space = vis.TextStim(win, text='Press space to close the experiment', pos=(0, -50), height=20)

# Draw final message
break_text.draw(); win.flip(); core.wait(5)  

# Draw both messages and wait for participant to press space
break_text.draw(); space.draw(); win.flip(); event.waitKeys(keyList=['space'])

        
# ---------------------------------------------------------------
# Save experimental data
# ---------------------------------------------------------------

# Save data as wide-format CSV
thisExp.saveAsWideText(file_name + '.csv', delim=',')

edf_filename = f"Data/RDM_reportz_eyetrack_{sub}.EDF"
tracker.setOfflineMode()
tracker.closeDataFile()
try:
    print("Receiving EDF from EyeLink...")
    tracker.receiveDataFile(edf_remote_name, edf_local_path)
    print(f"EDF saved to {edf_local_path}")
except RuntimeError as e:
    print("Error transferring EDF:", e)

tracker.close()


# ---------------------------------------------------------------
# End of experiment: close window and quit
# ---------------------------------------------------------------
win.close()
core.quit()


# ---------------------------------------------------------------
# Re-initialize ExperimentHandler (if starting new experiment/session)
# ---------------------------------------------------------------
file_name = "Data/RDM_reportz_sub%d" % sub
thisExp = data.ExperimentHandler(dataFileName=file_name,extraInfo=info)



