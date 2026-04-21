
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
import instructions_2 as ins
import EXP_functions as exp
from scipy.stats import norm
""" import pylink
import sys
import serial """ ###################################################################################################################################################################################################
#from pylink.eyelink import EyeLink
#from EyeLinkCoreGraphicsPsychoPy import EyeLinkCoreGraphicsPsychoPy


################################################################################################################################ 
# Set directory
################################################################################################################################ 
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

print(dname)
print(abspath)

################################################################################################################################ 
# Config
################################################################################################################################
# 
# ---------------------------------------------------------------
# Set up experiment parameters
# ---------------------------------------------------------------

pilot = 1                 # Pilot mode ON (1) or OFF (0)

# Staircase settings
SC_dotlife = 1            # Whether to run dotlife staircase (1 = yes, 0 = no)
SC_coherence = 1          # Whether to run coherence staircase (1 = yes, 0 = no)

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
resp_t = 0.5              # maximum time to respond
pres_time = 0.3        # presentation time during SC

# Accuracy values we want to present in the experiment
accVals = np.asarray([0.51, 0.6, 0.7, 0.8, 0.99])

# Number of repetitions for each coherence value in a sequence
repeats = [1, 2, 3, 2, 1]

# Maximum number of consecutive repetitions for up/down trials
max_cons = 4


# ---------------------------------------------------------------
# Trial characteristics
# ---------------------------------------------------------------

n_tutorial = 10 # Number of trials in second tutorial block

n_seq = 2                           # Number of times to present the sequence
len_seq = 2 * sum(repeats)          # Length of one sequence (sum of repeats * 2)

# Total number of trials per testing block
n_trials = sum(repeats) if pilot else n_seq * len_seq

# Total number of blocks
n_blocks = 9

# ---------------------------------------------------------------
# Hyperparameters for staircase dotlife
# ---------------------------------------------------------------

n_SC1 = 40           # Number of trials in dotlife staircase
coherence_d = 0.3    # Fixed coherence to find dotlife threshold performance

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

file_name = "Data_RDM2/RDM_metaWindow_sub%d" % sub  # Main data file path
thisExp = data.ExperimentHandler(dataFileName=file_name,   extraInfo=info)


# ---------------------------------------------------------------
# Setup file names and path for physio
# ---------------------------------------------------------------
""" sys.path.append("Eyelink")
edf_remote_name = f"RDM2_{sub}.EDF"  ###################################################################################################################################################################################################
edf_local_path = os.path.join("Data", f"RDM_metaWindow_eye{sub}.EDF") """


################################################################################################################################ 
# Psychopy objects
################################################################################################################################ 
# ---------------------------------------------------------------
# Set up Psychopy window
# ---------------------------------------------------------------
size = [1536, 960] if fullscreen else [600, 400]  # Window size depending on fullscreen mode
win = vis.Window(size=size,units='pix',color='grey',allowGUI=False,fullscr=fullscreen, waitBlanking=False)

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
choice_keys = ['1', '2', '3', '6', '7', '8', 'escape']  # for conf
resp_keys = ['3', '6', 'escape']    # for response

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

# Slider instructions
slider_instructions = vis.TextStim(win, text="How confident were you in your decision?", pos=(0, 75))


# ---------------------------------------------------------------
# Text for no-response feedback
# ---------------------------------------------------------------
no_response_text = vis.TextStim(win,text="No response, try to be faster next trial!",color='white',height=30)

# ---------------------------------------------------------------
# Create fixation cross and left right visuals
# ---------------------------------------------------------------
fixation = vis.ShapeStim(win=win,vertices=((0, -15), (0, 15), (0, 0), (-15, 0), (15, 0)),  # Cross arms
    lineWidth=5,closeShape=False,lineColor='white')

left = vis.TextStim(win,text= "L",pos=(-700, 0) ,height= 70,color='white',wrapWidth=1400,alignText='center')
right = vis.TextStim(win,text= "R",pos=(700, 0) ,height= 70,color='white',wrapWidth=1400,alignText='center')

# ---------------------------------------------------------------
# Space bar prompt for training / instructions
# ---------------------------------------------------------------
space = vis.TextStim(win, text='Press SPACE to continue', pos=(0, -300), height=30)


################################################################################################################################ 
# Set up Physio (copy paste from cobe lab set up)
################################################################################################################################ 
""" ###################################################################################################################################################################################################
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
 """

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
        mapping = {0: ('left', 180), 1: ('right', 0)}
        correct, direction = mapping[condition_direction[trial]]    

        # Prepare for trial
        resp = None          # Empty variable for response
        event.clearEvents()  # Clear any lingering key presses

        # Assign stimulus parameters
        DotMotion.coherence = coherence   # constant coherence
        DotMotion.dir = direction      # Set motion direction

        # ---------------------------------------------------------------
        # Present stimulus and wait for response
        # ---------------------------------------------------------------
        while not resp:
            fixation.draw()     # Draw fixation cross
            DotMotion.draw()    # Draw random dot stimulus
            win.flip()
            resp = event.getKeys(keyList=resp_keys)  # Collect response

        # ---------------------------------------------------------------
        # Calculate accuracy
        # ---------------------------------------------------------------
        correct_key = resp_keys[0] if correct == "left" else resp_keys[1]
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
        key_to_label = {resp_keys[0]: "up", resp_keys[1]: "down"}
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

    #index for the while loop
    trial = 0
    
    # Balanced left/right sequence for training
    sequence = n_SC1 // 2 * [0] + n_SC1 // 2 * [1]
    condition_direction = exp.generate_valid_sequence(sequence, max_cons)  # Avoid too many repeats

    # Initialize empty lists for accuracy and reaction times
    acc = [0] * n_SC1
    rt = [0] * n_SC1

    # Loop through all trials
    while trial < n_SC1:
        
        print("trial start")
        event.clearEvents()  # Clear keyboard events at start
        
        stim = SC.next()

        print("extracting SC value successful")
        
        # Map 0/1 to motion direction and angle
        mapping = {0: ('left', 180), 1: ('right', 0)}
        correct, direction = mapping[condition_direction[trial]]    
        

        # Assign stimulus parameters
        DotMotion.coherence = coherence_d   # constant coherence
        DotMotion.dotLife = stim          # dotLife set by staircase
        print(DotMotion.dotLife)
        DotMotion.dir = direction         # Set motion direction

        # ---------------------------------------------------------------
        # Present stimulus and wait for response
        # ---------------------------------------------------------------
        T_stimulus_start = clock.getTime()
        
        while (clock.getTime() - T_stimulus_start) < pres_time:
            fixation.draw()     # Draw fixation cross
            DotMotion.draw()    # Draw random dot stimulus
            win.flip()

        event.clearEvents()
        T_stimulus_end = clock.getTime()    
        while (clock.getTime() - T_stimulus_end) < resp_t:
            fixation.draw()
            left.draw()
            right.draw()
            win.flip()
            resp = event.getKeys(keyList=resp_keys)
            if resp:
                break
        
        # Timeout: skip trial if no response
        if not resp:
            print("No response within 700 ms, skipping trial")
            FB_text = "No response"
            FB_col = "white"
            RTdec = np.nan
            rt[trial] = RTdec
            ACC = 0 
            resp = None
        else: 
        # ---------------------------------------------------------------
        # Record reaction time
        # ---------------------------------------------------------------
            T_stimulus_stop = clock.getTime()
            RTdec = T_stimulus_stop - T_stimulus_end
            rt[trial] = RTdec
            print("response given")
            print("Reaction time is:", RTdec)

        # ---------------------------------------------------------------
        # Calculate accuracy
        # ---------------------------------------------------------------
        if resp:
            correct_key = resp_keys[0] if correct == "left" else resp_keys[1]
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
        key_to_label = {resp_keys[0]: "left", resp_keys[1]: "right"}
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

        trial +=1



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

    trial = 0

    # Loop through all trials
    while trial < n_SC2:

        # Get staircase value for this trial (coherence)
        stim = SC.next()
        print("Extracting SC value successful")

        # Map 0/1 to motion direction and angle
        mapping = {0: ('left', 180), 1: ('right', 0)}
        correct, direction = mapping[condition_direction[trial]]

        # Prepare trial
        event.clearEvents()  # Clear previous keyboard events

        # Assign stimulus parameters
        DotMotion.coherence = stim   # Coherence set by staircase
        DotMotion.dotLife = dotLife  # DotLife fixed
        DotMotion.dir = direction
        print(stim)                    # Print fcoherence 
        print("Presenting stimulus...")  # Print trial start

        # ---------------------------------------------------------------
        # Present stimulus and wait for response
        # ---------------------------------------------------------------
        T_stimulus_start = clock.getTime()
        
        while (clock.getTime() - T_stimulus_start) < pres_time:
            fixation.draw()     # Draw fixation cross
            DotMotion.draw()    # Draw random dot stimulus
            win.flip()

        event.clearEvents()
        T_stimulus_end = clock.getTime()    
        while (clock.getTime() - T_stimulus_end) < resp_t:
            fixation.draw()
            left.draw()
            right.draw()
            win.flip()
            resp = event.getKeys(keyList=resp_keys)
            if resp:
                break

        if not resp:
            print("No response within 700 ms, skipping trial")
            FB_text = "No response"
            FB_col = "white"
            RTdec = np.nan
            rt[trial] = RTdec
            ACC = 0 
            resp = None
        else: 
        # ---------------------------------------------------------------
        # Record reaction time
        # ---------------------------------------------------------------
            T_stimulus_stop = clock.getTime()
            RTdec = T_stimulus_stop - T_stimulus_end
            rt[trial] = RTdec
            print("response given")
            print("Reaction time is:", RTdec)


        # ---------------------------------------------------------------
        # Calculate accuracy
        # ---------------------------------------------------------------
            correct_key = resp_keys[0] if correct == "left" else resp_keys[1]
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
        key_to_label = {resp_keys[0]: "left", resp_keys[1]: "right"}
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

        trial += 1

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
    ins.Main1(win); core.wait(1); event.waitKeys(keyList=['lctrl'])
    ins.Main2(win); core.wait(1); event.waitKeys(keyList=['space'])
    ins.Main3(win); core.wait(1); event.waitKeys(keyList=['space'])
    ins.Main4(win); core.wait(1); event.waitKeys(keyList=['space'])



thisExp.saveAsWideText(file_name + '.csv', delim=',')
win.close()
core.quit()