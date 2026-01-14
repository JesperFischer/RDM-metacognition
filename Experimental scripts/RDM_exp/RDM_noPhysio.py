# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 09:47:57 2025
@author: Siebe Everaerts

Code for presenting random dot motion stimuli in 2 directions, recording responses
and giving feedback on the decision. Saves reaction time, response, accuracy, and
participant characteristics.

Use by pressing the "q" (or "a") key when dots are moving left and the "e" key when dots are moving right

"""

# Set up experiment -----------------------------------------------------------
# Import modules
import random
import numpy as np
from psychopy import visual as vis
from psychopy import event, core, data
from psychopy.hardware import keyboard
from psychopy.data import QuestPlusHandler
from scipy.stats import truncnorm
import questplus as qp
import concurrent.futures
import math
import itertools
import os
#import pylink
import sys
import instructions_1 as ins

#sys.path.append("Eyelink")
#ser = serial.Serial('COM7', 115200, timeout=1)


################################################################################################################################ 
# Set directory
################################################################################################################################ 
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)



################################################################################################################################ 
# Experiment variables
################################################################################################################################ 
# general
keyboard_type = "qwerty" # azerty
pilot = 1
SC_dotlife = 1
training = 0
training_2 = 1
instructions = 1
fullscreen = 1 #Developper mode
Part = 0 # 0 = Part 1 and 1 = Part 2


# training
if not pilot: 
    n_training = 20 # number of training trials per block
    n_training_2 = 10 
else:
    n_training = 5
    n_training_2 = 4
  
per_correct = 0 # percentage correct when starting --> for training if percetage is lower start another block 
des_per_cor = 0.6 # desired percentage correct
mean_rt = 5 # mean reaction time before starting (seconds)
des_mean_rt = 2 # desired mean reaction time --> for training if percetage is lower start another block 
coherence = .80 # coherence for first training block 
coherence_hard = .40 # coherence for second training block (staircase sets it for main exp)
max_dur_conf = 3 # maximum duration for confidence scale 

# main blocks
if not pilot:
    n_trials = 46 # Number of trials per testing block Note: this should be an even number
    n_blocks = 8 # Number of testing blocks 
else:
    n_trials = 10
    n_blocks = 4

block = 0 # starting block number (for training/staircase)

# stimulus
dotLife = 5 

# waiting times (in seconds)
if not pilot:
    ins_wait = 1; break_wait = 5
elif pilot:
    ins_wait = 0; break_wait = 0

# Determine when other scale is presented
if not pilot:
    scale_Nblock = 4
elif pilot:
    scale_Nblock = 1


################################################################################################################################ 
# Setting up a data file
################################################################################################################################ 
#if pilot:
#    sub = 0; age = 30; gender = 'Man'; handedness = 'Right'
#    info = {"sub": sub, "age": age, "gender": gender, "handedness": handedness} #creating a dictionnary
#else:
#    sub = int(input("Subject number: "))
#    age = int(input("Age: "))
#    gender = input("Gender (Woman/Man/X): ")
#    handedness = input("Handedness (Left/Right): ")
#    info = {"subject": sub, "age": age, "gender": gender,"handedness": handedness}
  

sub = 0; age = 30; gender = 'Man'; handedness = 'Right'
info = {"sub": sub, "age": age, "gender": gender, "handedness": handedness} #creating a dictionnary


file_name = "Data/RDM_reportz_sub%d" % sub
thisExp = data.ExperimentHandler(dataFileName=file_name, extraInfo=info)  # saving extra info along with the main experimental data

## file save for eyelink 1000
# remote filename (on the EyeLink) must be a simple filename (no path)
#edf_remote_name = f"sqq_{sub}.EDF"
# local path where we'll save the received EDF
#edf_local_path = os.path.join("Data", f"RDM_reportz_eyetrack_{sub}.EDF")


################################################################################################################################ 
# Psychopy objects
################################################################################################################################ 

if fullscreen:
    win = vis.Window(size=[1536,960], units = 'pix', color='grey', allowGUI=False, fullscr=True) #set fullscr to True for the experiment
else:
    win = vis.Window(size = [600,400], units = 'pix', color='grey', allowGUI=False, fullscr=False)


width = win.size[0]
height = win.size[1]

#Create a keybard object and stepsize for slider
kb = keyboard.Keyboard()
step_size = 0.02 

# Clock
clock = core.Clock()

# Define response keys based on keyboard type
if keyboard_type == "qwerty":
    choice_keys = ['w', 's', 'escape']  # up, down, escape
elif keyboard_type == "azerty":
     choice_keys = ['z', 's', 'escape']  # up, down, escape
else:
     raise TypeError('Unknown keyboard name')

# Creating DotMotion stimulus
DotMotion = vis.DotStim(win, units='pix', nDots= 120, fieldSize = 300, fieldShape='circle', dotSize=7  , dotLife=15, speed=1.5, color='white', 
                        signalDots='same', noiseDots='walk') #https://www.psychopy.org/api/visual/dotstim.html 

# Creating a slider to rate confidence or clarity 
slider = vis.Slider(win, name='slider', size=(400,20), pos = (0,0), units = 'pix',
                          ticks=(0,1), granularity = 0.01,
                          style=['rating'], color='white', font='HelveticaBold', flip=False)
slider.marker.color = "white"
slider.marker.size = 20                       
slider_label_wrong = vis.TextStim(win, text= "definitely wrong", pos=(-200, 30)) 
slider_label_right = vis.TextStim(win, text= "definitely right", pos=(200, 30))                     
slider_label_Nclear = vis.TextStim(win, text= "not clear at all", pos=(-200, 30)) 
slider_label_clear = vis.TextStim(win, text= "very clear", pos=(200, 30)) 
slider_instructions = vis.TextStim(win, text = "How confident were you in your decision?", pos=(0,75))
slider_instructions_dir = vis.TextStim(win, text = "How clear was the DIRECTION of the dots?", pos=(0,75))
mouse = event.Mouse(win=win) #Create a mouse object for the slider

#Text for confidence no response
no_response_text = vis.TextStim(win, text="No response, try to be faster next trial!",
            color='white', height=30)

# creating a fixation cross
fixation = vis.ShapeStim(
    win=win,
    vertices=((0, -15), (0, 15), (0,0), (-15,0), (15,0)),  # vertical and horizontal lines
    lineWidth=5,
    closeShape=False,
    lineColor='white'
)

# space for training
space = vis.TextStim(win, text='Press SPACE to continue', pos=(0, -300), height=30)

# === Launch iohub with EyeLink ===
#tracker = EyeLink("100.1.1.1")
#tracker.openDataFile(edf_remote_name)

# === EyeLink calibration graphics ===
#genv = EyeLinkCoreGraphicsPsychoPy(tracker, win)
#genv.setCalibrationColors((-1, -1, -1), win.color)
#genv.setTargetType('picture')
#genv.setPictureTarget(os.path.join(dname, 'Eyelink_cal', 'fixTarget.bmp'))
#pylink.openGraphicsEx(genv)

# === EyeLink screen and calibration settings ===
#tracker.sendCommand(f"screen_pixel_coords = 0 0 {width - 1} {height - 1}")
#tracker.sendMessage(f"DISPLAY_COORDS = 0 0 {width - 1} {height - 1}")
#tracker.sendCommand('enable_automatic_calibration=YES')
#tracker.sendCommand('automatic_calibration_pacing=500')
           
##################################################################
######## Eyetracker  Calibration #################################
##################################################################
               
# === Calibrate the tracker (optional but recommended) ===
#tracker.doTrackerSetup()






# reading instructions slides
#intro = vis.ImageStim(win, image=dname+"\Intro.jpg", units = 'pix', size = [width,height])
#main_1 = vis.ImageStim(win, image=dname+"\Main1.jpg", units = 'pix', size = [width, height]) 
#main_2 = vis.ImageStim(win, image=dname+"\Main2.jpg", units = 'pix', size = [width, height]) 
#main_3 = vis.ImageStim(win, image=dname+"\Main3.jpg", units = 'pix', size = [width, height])
#main_4 = vis.ImageStim(win, image=dname+"\Main4.jpg", units = 'pix', size = [width, height])
#main_5 = vis.ImageStim(win, image=dname+"\Main5.jpg", units = 'pix', size = [width, height]) 
#main_6 = vis.ImageStim(win, image=dname+"\Main6.jpg", units = 'pix', size = [width, height])     

################################################################################################################################ 
# Functions
################################################################################################################################ 
def get_stim(SC):
    return SC.next_stim

## Function to check if the mouse is hovering over the slider bar area
def move_slider(left, right, up, slider, SR, step_size):
    slider_pos = slider.markerPos
    if slider_pos is None:
        slider_pos = 0.5
    if left:
        slider_pos = max(0, slider.markerPos - step_size)  
    if right:
        slider_pos = min(1, slider.markerPos + step_size)
    if up:
        slider_pos = slider.markerPos
        SR = slider_pos

    return slider_pos, SR 

# Get training state
def get_state(per_correct, mean_rt, des_per_cor, des_mean_rt, TrialType):
    result = {}

    if per_correct >= des_per_cor and mean_rt <= des_mean_rt and TrialType == "Training - easy":
        result['feedback'] = "Great job! We're going to make it a little more difficult now."
        result['state'] = 0

    elif per_correct >= des_per_cor and mean_rt <= des_mean_rt and TrialType == "Training - hard":
        result['feedback'] = None
        result['state'] = 1

    elif per_correct <= des_per_cor or mean_rt > des_mean_rt:
        result['percentage'] = f'You got {per_correct*100:.0f}% correct.'
        result['time'] = f'Your average reaction time was {mean_rt:.2f} seconds'

        if per_correct < des_per_cor and mean_rt > des_mean_rt:
            result['feedback'] = "Try to be faster and more accurate in the next block!"
            result['state'] = 2
        elif per_correct < des_per_cor:
            result['feedback'] = "Try to be more accurate in the next block!"
            result['state'] = 3
        elif mean_rt > des_mean_rt:
            result['feedback'] = "Try to be faster in the next block!"
            result['state'] = 4
    else:
        result['feedback'] = None
        result['state'] = None

    return result

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

################################################################################################################################ 
# Training 1
################################################################################################################################ 
win.mouseVisible = False

# welcome text
ins.Intro(win);core.wait(ins_wait); event.waitKeys(keyList=['space'])

# training blocks
if training:
    TrialType = "Training - easy"
    while (per_correct <= des_per_cor or mean_rt > des_mean_rt):
        print("Conditions not met, starting another practice block")

        # training: 50% left and right
        condition_direction = np.repeat(range(2),[math.floor(n_training*0.5), math.ceil(n_training*0.5)]); random.shuffle(condition_direction) #Creates an equal amount of left/right trials

        #Empty lists of accuracy and reaction times for training data  
        acc = [0] * n_training 
        rt = [0] * n_training
        for trial in range(n_training):
            # Stimulus direction
            mapping = {0: ('up', 90), 1: ('down', 270)}
            correct, direction = mapping[condition_direction[trial]]    

            # draw stimulus
            resp = None #empty list for response
            event.clearEvents() 
            DotMotion.coherence = coherence
            DotMotion.dotLife = dotLife
            DotMotion.dir = direction

            # save start time of the stimulus    
            T_stimulus_start = clock.getTime()
            while not resp:
                fixation.draw()
                DotMotion.draw()
                win.flip()
                resp = event.getKeys(keyList=choice_keys)
                if clock.getTime() - T_stimulus_start >= des_mean_rt:
                    print("No response within 2 s, skipping trial")
                    FB_text = "No response"
                    FB_col = "white"
                    break
                       
            if resp:
                T_stimulus_stop = clock.getTime()
                RTdec = T_stimulus_stop - T_stimulus_start
                rt[trial] = RTdec
                print("Reaction time is:", RTdec)
            else:
                RTdec = np.nan
                rt[trial] = RTdec
        
            
            # get accuracy
            correct_key = choice_keys[0] if correct == "up" else choice_keys[1]
            if resp:
                is_correct = (resp[0] == correct_key)
                ACC = int(is_correct)
                acc[trial] = ACC
    
                if is_correct:
                    print("Decision was correct")
                    FB_text = "Correct!"
                else:
                    print("Decision was incorrect")
                    FB_text = "Wrong"
        
            else:
                 ACC = 0
            
            # allow escape to exit experiment
            if resp == ['escape']:
                print('Participant pressed escape')
                thisExp.saveAsWideText(file_name + '.csv', delim=',') 
                win.close()
                core.quit()
            
            # Give feedback
            feedback = vis.TextStim(win, text = FB_text, color = "white", height=40)
            feedback.draw()
            win.flip()

            core.wait(0.5 if resp else 1)

            key_to_label = {choice_keys[0]: "up", choice_keys[1]: "down"}

            if resp:
                resp = key_to_label[resp[0]]
                
            thisExp.addData("block", block)
            thisExp.addData("Trialtype", TrialType)
            thisExp.addData("withinblocktrial", trial)
            thisExp.addData("RTdec", RTdec)
            thisExp.addData("resp", resp)
            thisExp.addData("cor", ACC)
            thisExp.addData("dots direction", direction)
            thisExp.addData("cor_resp", correct)
            thisExp.addData("coherence", DotMotion.coherence)
            thisExp.nextEntry()
            
        per_correct = sum(acc) / len(acc)
        mean_rt = np.nanmean(rt)
        print("Percent correct is:", per_correct)
        print("Average reaction time is:", mean_rt)

        state = get_state(per_correct, mean_rt, des_per_cor, des_mean_rt, TrialType)

        if state["state"] == 0:
            print("Conditions met, moving on to harder training block")
            
            feedback_text = vis.TextStim(win, text = state["feedback"], pos = (0, 200), height = 40, wrapWidth = 1200)
            break_text = vis.TextStim(win, text = "Take a short break before we continue with the next block.", pos = (0,0), height = 40, wrapWidth = 1200)

            feedback_text.draw(); break_text.draw(); win.flip()
            core.wait(break_wait); feedback_text.draw(); break_text.draw(); space.draw(); win.flip(); event.waitKeys(keyList=['space']) #Use corewait for intertrial interval

            TrialType = "Training - hard"
            coherence = coherence_hard
            per_correct = 0; mean_rt = 2; # set back to default for new run
             
        elif state["state"] == 1:
            print("Conditions met, moving on to real experiment")
            break
        
        elif state["state"] > 1: 
               
            # offer a break
            points_text = vis.TextStim(win, text = state["percentage"], pos=(0, 200),height = 35, wrapWidth = 1200)
            speed_text = vis.TextStim(win, text = state["time"], pos = (0,150), height = 35, wrapWidth = 1200)
        
            if state["state"] == 2:
                feedback_text = vis.TextStim(win, text = state["feedback"], height = 35, wrapWidth = 1200, pos = (0,50))
            elif state["state"] == 3:
                feedback_text = vis.TextStim(win, text = state["feedback"], height = 35, wrapWidth = 1200, pos = (0,50))
            elif state["state"] == 4:
                feedback_text = vis.TextStim(win, text = state["feedback"], height = 35, wrapWidth = 1200, pos = (0,50))
        
            break_text = vis.TextStim(win, text = "Take a short break before we continue with the next block.", pos = (0,-50), height = 35, wrapWidth = 1200)
            space = vis.TextStim(win, text='Press space to continue', pos=(0, -250), height=30)
            points_text.draw(); speed_text.draw(); feedback_text.draw(); break_text.draw(); win.flip()
            core.wait(break_wait); points_text.draw(); speed_text.draw(); feedback_text.draw(); break_text.draw(); space.draw(); win.flip(); event.waitKeys(keyList=['space']) 

################################################################################################################################ 
# Staircase Dotlife
################################################################################################################################         
TrialType = "SC dotlife" # Trialtype
block += 1
n_SC1 = 40   # number of trials in staircase
coherence = 0.3 # Set fixed value to find dotlife (task difficulty) where threshold performance should be at 0.3 given a dotLife

#alpha = 5
#beta = 1.3
#def weibull_cdf(x, alpha, beta, gamma =0.5, lapse =0):
    #return 1 - lapse - (1 - gamma - lapse) * (np.exp(-(x / alpha)**beta))
 


# STAIRCASE (https://questplus.readthedocs.io/en/latest/qp.html)
SC = QuestPlusHandler(nTrials = n_SC1, intensityVals = list(range(0, 101, 5)), 
                      thresholdVals = np.linspace(0,100,50), slopeVals=np.linspace(1,10,10),lowerAsymptoteVals = 0.5, lapseRateVals=0,
                      responseVals = [1,0], prior = {"threshold": np.ones(50)/50, "slope": np.ones(10)/10}, psychometricFunc = "weibull", startIntensity = 30,
                      stimScale="linear", stimSelectionMethod="minEntropy", paramEstimationMethod = "mean")

if SC_dotlife:
        # training: 50% left and right
        condition_direction = np.repeat(range(2),[math.floor(n_SC1*0.5), math.ceil(n_SC1*0.5)]); random.shuffle(condition_direction) #Creates an equal amount of left/right trials

        #Empty lists of accuracy and reaction times for training data  
        acc = [0] * n_SC1 
        rt = [0] * n_SC1
        for trial in range(n_SC1):
            stim = SC.next()
            # Stimulus direction
            mapping = {0: ('up', 90), 1: ('down', 270)}
            correct, direction = mapping[condition_direction[trial]]    

            # draw stimulus
            resp = None #empty list for response
            event.clearEvents() 
            DotMotion.coherence = coherence
            DotMotion.dotLife = stim
            print(DotMotion.dotLife)
            DotMotion.dir = direction

            # save start time of the stimulus    
            T_stimulus_start = clock.getTime()
            while not resp:
                fixation.draw()
                DotMotion.draw()
                win.flip()
                resp = event.getKeys(keyList = choice_keys)
                #resp = "up"
                if clock.getTime() - T_stimulus_start >= des_mean_rt:
                    print("No response within 2 s, skipping trial")
                    FB_text = "No response"
                    FB_col = "white"
                    break
                       
            if resp:
                T_stimulus_stop = clock.getTime()
                RTdec = T_stimulus_stop - T_stimulus_start
                rt[trial] = RTdec
                print("Reaction time is:", RTdec)
            else:
                RTdec = np.nan
                rt[trial] = RTdec
        
            
            #get accuracy
            correct_key = choice_keys[0] if correct == "up" else choice_keys[1]
            if resp:
                is_correct = (resp[0] == correct_key)
                ACC = int(is_correct)
                acc[trial] = ACC
            #prob = weibull_cdf(stim, alpha, beta)
            #print(prob)
            #ACC = np.random.binomial(n=1, p=prob)
                if is_correct:
                    print("Decision was correct")
                    FB_text = "Correct!"
                else:
                    print("Decision was incorrect")
                    FB_text = "Wrong"
        
            else:
                 ACC = 0
            
            # allow escape to exit experiment
            if resp == ['escape']:
                print('Participant pressed escape')
                thisExp.saveAsWideText(file_name + '.csv', delim=',') 
                win.close()
                core.quit()

            # Update staircase
            print("trial = OK")
            SC.addResponse(ACC)
            print("update = OK")
        
            # Give feedback
            feedback = vis.TextStim(win, text = FB_text, color = "white", height=40)
            feedback.draw()
            win.flip()

            core.wait(0.5 if resp else 1)

            key_to_label = {choice_keys[0]: "up", choice_keys[1]: "down"}

            if resp:
                resp = key_to_label[resp[0]]
                
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

print(SC.paramEstimate["threshold"])
dotLife = SC.paramEstimate["threshold"]

################################################################################################################################ 
# Instructions
################################################################################################################################ 

if instructions:
    ins.Main1(win);core.wait(ins_wait); event.waitKeys(keyList=['space'])  
    ins.Main2(win);core.wait(ins_wait); event.waitKeys(keyList=['space']) 
    ins.Main3(win);core.wait(ins_wait); event.waitKeys(keyList=['space']) 
    ins.Main4(win);core.wait(ins_wait); event.waitKeys(keyList=['space'])  
    ins.Main5(win);core.wait(ins_wait); event.waitKeys(keyList=['space']) 

################################################################################################################################ 
# Training 2
################################################################################################################################
block += 1
TrialType = "Training - scale"
if training_2:
    #draw intensities
    coherence = np.linspace(0.1, 1.0, n_training_2);random.shuffle(coherence)

    #draw directions
    condition_direction = np.repeat(range(2),[math.floor(n_training_2*0.5), math.ceil(n_training_2*0.5)]); random.shuffle(condition_direction)

    acc = [0] * n_training_2
    rt = [0] * n_training_2
    for trial in range(n_training_2):
        # Stimulus direction
        mapping = {0: ('up', 90), 1: ('down', 270)}
        correct, direction = mapping[condition_direction[trial]]  
                
        # draw stimulus
        resp = None #empty list for response
        event.clearEvents() 
        DotMotion.coherence = coherence[trial]
        DotMotion.dotLife = dotLife
        DotMotion.dir = direction

        #stimulus loop   
        T_stimulus_start = clock.getTime()
        while not resp:
            fixation.draw()
            DotMotion.draw()
            win.flip()
            resp = event.getKeys(keyList=choice_keys)
            if clock.getTime() - T_stimulus_start >= des_mean_rt:
                print("No response within 1.5 s, skipping trial")
                miss_text = vis.TextStim(win, text = "No response, try to be faster next trial!", height = 30)
                miss_text.draw(); win.flip()
                core.wait(1)
                break
                    
        if resp:
            T_stimulus_stop = clock.getTime()
            RTdec = T_stimulus_stop - T_stimulus_start
            rt[trial] = RTdec
            print("Reaction time is:", RTdec)
        else:
            RTdec = np.nan
            rt[trial] = RTdec
        
        # get accuracy
        correct_key = choice_keys[0] if correct == "up" else choice_keys[1]
        if resp:
            is_correct = (resp[0] == correct_key)
            ACC = int(is_correct)
            acc[trial] = ACC
        else:
                ACC = 0

        # allow escape to exit experiment
        if resp == ['escape']:
            print('Participant pressed escape')
            thisExp.saveAsWideText(file_name + '.csv', delim=',') 

            win.close()
            core.quit()
        
        fixation.draw(); win.flip(); core.wait(1)
        kb.clearEvents()


        if resp:
            T_rating_start = clock.getTime()
            slider.reset()
            slider.draw()
            slider_instructions.draw(); slider_label_wrong.draw(); slider_label_right.draw()    
            win.flip()

            SR  = None
            held_keys =[]
            while SR is None: 
                # check if participant is 
                elapsed_time = clock.getTime() - T_rating_start
                if elapsed_time >= max_dur_conf:
                    no_response_text.draw(); win.flip()
                    core.wait(1)
                    SR = None  # make sure response is recorded as missing
                    RTrating = None
                    break

                # Get Key presses
                left, right, up, escape = kb.getState(['left', 'right', 'up', 'escape'])
                
                # Check if escape
                if escape:  
                    print('Participant pressed escape')
                    thisExp.saveAsWideText(file_name + '.csv', delim=',') 
                    win.close()  
                    core.quit()

                # Get new slider value and check if new confidence
                slider_pract_value, SR = move_slider(left, right, up, slider, SR, step_size)
                slider.markerPos = slider_pract_value  # Update slider marker position
                slider.draw()
                slider_instructions.draw(); slider_label_wrong.draw(); slider_label_right.draw()    

                win.flip()

                # Check if the mouse has been clicked to submit the answer
            print("Reported confidence = ", SR)
            T_rating_stop = clock.getTime()
            RTrating = T_rating_stop - T_rating_start
            
        else:
            RTrating = None
            SR = None

        fixation.draw(); win.flip(); core.wait(1)

        key_to_label = {choice_keys[0]: "up", choice_keys[1]: "down"}

        if resp:
            resp = key_to_label[resp[0]]
        else:
            scale = None

        # Save trial
        thisExp.addData("block", block)
        thisExp.addData("Trialtype", TrialType)
        thisExp.addData("withinblocktrial" , trial)
        thisExp.addData("RTdec", RTdec)
        thisExp.addData("resp", resp)
        thisExp.addData("cor", ACC)
        thisExp.addData("dots direction", direction)
        thisExp.addData("cor_resp", correct)
        thisExp.addData("SR_conf", SR)
        thisExp.addData("RTrating", RTrating)
        thisExp.addData("coherence", coherence[trial])
        thisExp.nextEntry()

ins.Main6(win);core.wait(ins_wait); event.waitKeys(keyList=['space'])

################################################################################################################################ 
# Variables for main experiment
################################################################################################################################
TrialType = "Main" # Trialtype
inter_t_mean = 2.5; inter_t_sd = 1; lower = 1; upper = 4 # variables inter-trial-interval
a = standard(lower,inter_t_mean,inter_t_sd); b = standard(upper, inter_t_mean,inter_t_sd) # standardized for truncated normal
means = [2.375, 4.125]; sds = [0.05, 0.5] #variables distributions Part 2 
pairs = list(itertools.product(means, sds)); repeat = n_blocks // len(pairs); conditions = repeat * pairs; np.random.shuffle(conditions) #conditions for Part 2

# STAIRCASE (https://questplus.readthedocs.io/en/latest/qp.html)
SC = qp.QuestPlus(stim_domain= {"intensity": np.linspace(0.01,1,50)}, 
                 func="weibull",
                 stim_scale="log10",
                 param_domain= {"threshold": np.linspace(0.01,1,50), "slope": np.linspace(1,10,50), "lower_asymptote": 0.5, "lapse_rate": 0.05},
                 prior = {"threshold": np.ones(50)/50, "slope": np.ones(50)/50},
                 outcome_domain={"response": [1, 0]},
                 stim_selection_method="min_n_entropy",
                 stim_selection_options = {"n": 1, "max_consecutive_reps": 20},
                 param_estimation_method= "mean")

# Equal # trials left and right for first block 
condition_direction = np.repeat(range(2),[math.floor(n_trials*0.5), math.ceil(n_trials*0.5)]); random.shuffle(condition_direction)

# determine waiting times between trials for first block
inter_trial = truncnorm.rvs(a, b, loc=inter_t_mean, scale=inter_t_sd, size= n_trials) #Can change mean according to pilots

#waiting times confidence interval for first block
if not Part:
    manipulation = np.random.uniform(low = des_mean_rt, high = 5, size = n_trials) #Based on Bradley et al. (2012) "Orienting and Emotional Perception: Facilitation, Attenuation, and Interference"
else:
    current_mean = conditions[0][0]; current_sd = conditions[0][1] #current variables for normal distr
    a2 = standard(1.5,current_mean, current_sd); b2 = standard(5, current_mean, current_sd) #normalize
    manipulation = truncnorm.rvs(a2, b2, loc= current_mean, scale=current_sd, size= n_trials) 

#Initiate vectors for first block 
acc = [0] * n_trials
rt = [0] * n_trials
#Trial number and block
trialN = 0
blockN = 0

################################################################################################################################ 
# Main experiment
################################################################################################################################ 
for eachTrial in range(n_trials*n_blocks):

    #tracker.setOfflineMode()

    # Stimulus direction
    mapping = {0: ('up', 90), 1: ('down', 270)}
    correct, direction = mapping[condition_direction[trialN]]
            
    # save start time of the stimulus    
    T_stimulus_start = clock.getTime()
    
    # draw stimulus
    with concurrent.futures.ThreadPoolExecutor() as executor:
        future = executor.submit(get_stim, SC)
        try:
            coherence = future.result(timeout=2)  # seconds
        except concurrent.futures.TimeoutError:
            print(f"Problem on trial: {trialN}, block: {blockN + 2}. Aborting experiment...")
            thisExp.saveAsWideText(file_name + '.csv', delim=',')
            ## save eyelink EDF from tracker to local Data folder
            #tracker.setOfflineMode()
            #tracker.closeDataFile()
            try:
                print("Receiving EDF from EyeLink...")
                #tracker.receiveDataFile(edf_remote_name, edf_local_path)
                #print(f"EDF saved to {edf_local_path}")
            except RuntimeError as e:
                print("Error transferring EDF:", e)

            #tracker.close()
            win.close()
            core.quit()
    
    print(coherence["intensity"])
    resp = None
    event.clearEvents() 
    DotMotion.coherence = coherence["intensity"]
    DotMotion.dotLife = dotLife
    DotMotion.dir = direction

    ######## start the eyetracker:
    #tracker.startRecording(1, 1, 1, 1)
    #tracker.sendMessage(f"start_trialID_{trialN}_Block_{blockN + 2}")
    ### biopack
    #ser.write(str.encode('01'))
    #core.wait(0.1)
    #ser.write(str.encode('00')) # turn off all 8

    #tracker.sendMessage("start_stimulus")  # Optional: timestamp visual onset
    while not resp:
        #win.color = 'white' 
        fixation.draw()
        DotMotion.draw()
        win.flip()
        resp = event.getKeys(keyList=choice_keys)
        if clock.getTime() - T_stimulus_start >= des_mean_rt:
            print("No response within 1.5 s, skipping trial")  
            miss_text = vis.TextStim(win, text = "No response, try to be faster next trial!", height = 30)
            miss_text.draw(); win.flip()
            core.wait(1)
            break
    #tracker.sendMessage("end_stimulus")  # Optional: timestamp visual offset  
    # get reaction time
    if resp:
        T_stimulus_stop = clock.getTime()
        RTdec = T_stimulus_stop - T_stimulus_start
        rt[trialN] = RTdec
        print("Reaction time is:", RTdec)
    else:
        rt[trialN] = np.nan
        
    # get accuracy
    correct_key = choice_keys[0] if correct == "up" else choice_keys[1]
    if resp:
        is_correct = (resp[0] == correct_key)
        ACC = int(is_correct)
        acc[trialN] = ACC

    else:
            ACC = 0
            
    # allow escape to exit experiment
    if resp == ['escape']:
        print('Participant pressed escape')
        thisExp.saveAsWideText(file_name + '.csv', delim=',') 

        ## save eyelink EDF from tracker to local Data folder
        #tracker.setOfflineMode()
        #tracker.closeDataFile()
        try:
            print("Receiving EDF from EyeLink...")
            #tracker.receiveDataFile(edf_remote_name, edf_local_path)
            #print(f"EDF saved to {edf_local_path}")
        except RuntimeError as e:
            print("Error transferring EDF:", e)

        #tracker.close()

        win.close()
        core.quit()
    
    #fixation cross
    #win.color = 'black'
    fixation.draw() ; win.flip()
    kb.clearEvents()

    #Waiting time for cofidence ratings
    if resp:
        interval = manipulation[trialN]
        core.wait(interval - rt[trialN]- 0.05)

        #ser.write(str.encode('01'))
        #core.wait(0.05)
        #ser.write(str.encode('00')) # turn off all 8


        T_rating_start = clock.getTime()
        slider.reset()
        start_conf = np.random.uniform(low=0.25, high=0.75)
        slider.setMarkerPos(start_conf)
        #win.color = 'white'
        slider.draw()
        if blockN == scale_Nblock:
            slider_instructions_dir.draw(); slider_label_Nclear.draw(); slider_label_clear.draw()
        else:
            slider_instructions.draw(); slider_label_wrong.draw(); slider_label_right.draw()
        win.flip()


        SR = None
        #tracker.sendMessage("start_confidence")
        while SR is None: 
            # check if participant is 
            elapsed_time = clock.getTime() - T_rating_start
            if elapsed_time >= max_dur_conf:
                no_response_text.draw(); win.flip(); core.wait(1)
                SR = None  # make sure response is recorded as missing
                RTrating = None
                break
# Get Key presses 
            left, right, up, escape = kb.getState(['left', 'right', 'up', 'escape'])

            # Check if escape
            if escape:  
                print('Participant pressed escape')
                thisExp.saveAsWideText(file_name + '.csv', delim=',') 
                                ## save eyelink EDF from tracker to local Data folder
                #tracker.setOfflineMode()
                #tracker.closeDataFile()
                try:
                    print("Receiving EDF from EyeLink...")
                    #tracker.receiveDataFile(edf_remote_name, edf_local_path)
                    #print(f"EDF saved to {edf_local_path}")
                except RuntimeError as e:
                    print("Error transferring EDF:", e)
                #tracker.close()
                win.close()  
                core.quit()

            # Get new slider value and check if new confidence
            slider_pract_value, SR = move_slider(left, right, up, slider, SR, step_size)
            slider.markerPos = slider_pract_value  # Update slider marker position
        
            # Redraw the slider and instructions
            #win.color = 'white'
            slider.draw()
        
            if blockN == scale_Nblock:
                slider_instructions_dir.draw(); slider_label_Nclear.draw(); slider_label_clear.draw()
                scale = 'control'
            else:
                slider_instructions.draw(); slider_label_wrong.draw(); slider_label_right.draw() 
                scale = 'conf'
            win.flip()
            
            T_rating_stop = clock.getTime()
            RTrating = T_rating_stop - T_rating_start
        
    else:
        RTrating = None
        SR = None
        interval = None
    
    core.wait(0.1)
    #tracker.sendMessage(f"end_trialID_{trialN}_Block_{blockN + 2}")
    #ser.write(str.encode('01'))
    #core.wait(0.015)
    #ser.write(str.encode('00'))
    #tracker.stopRecording()


    #Add response to staircase and proceed to next value
    print("trial = OK")
    SC.update(stim= coherence, outcome= {"response":ACC})
    print("update = OK")
        

    # Blank screen drawn from a truncated normal distribution
    #win.color = 'black'
    fixation.draw(); win.flip(); core.wait(inter_trial[trialN]) # Change to waiting time drawn from a distribution 

    key_to_label = {choice_keys[0]: "up", choice_keys[1]: "down"}
    if resp:
        resp = key_to_label[resp[0]]
    else:
        scale = None
        
    # Save trial
    thisExp.addData("block", blockN + 2)
    thisExp.addData("Trialtype", TrialType)
    thisExp.addData("withinblocktrial" , trialN)
    thisExp.addData("RTdec", RTdec)
    thisExp.addData("resp", resp)
    thisExp.addData("cor", ACC)
    thisExp.addData("dots direction", direction)
    thisExp.addData("cor_resp", correct)
    thisExp.addData("interval", interval)
    thisExp.addData("interTrial.interval", inter_trial[trialN])
    if Part:
        thisExp.addData("Mean", conditions[blockN][0])
        thisExp.addData("Standard deviation", conditions[blockN][1])
    thisExp.addData("scale", scale)
    thisExp.addData("start_conf", start_conf)
    thisExp.addData("SR_conf", SR)
    thisExp.addData("RTrating", RTrating)
    thisExp.addData("coherence", coherence["intensity"])
    thisExp.nextEntry()

    #Update trialN
    trialN += 1

    if blockN < n_blocks - 1:
    #Check if next block     
        if trialN == n_trials: 
            blockN += 1
            trialN = 0

            #Update variables for next block
            condition_direction = np.repeat(range(2),[math.floor(n_trials*0.5), math.ceil(n_trials*0.5)]); random.shuffle(condition_direction)

            # determine waiting times between trials and waiting times confidence interval
            inter_trial = truncnorm.rvs(a, b, loc=inter_t_mean, scale=inter_t_sd, size= n_trials) #Can change mean according to pilots

            if not Part:
                manipulation = np.random.uniform(low = des_mean_rt, high = 5, size = n_trials) #Based on Bradley et al. (2012) "Orienting and Emotional Perception: Facilitation, Attenuation, and Interference"
            else:
                current_mean = conditions[blockN][0]; current_sd = conditions[blockN][1] # parameters for distribution
                a2 = standard(1.5,current_mean, current_sd); b2 = standard(5, current_mean, current_sd) # normalized parameters
                manipulation = truncnorm.rvs(a2, b2, loc= current_mean, scale=current_sd, size= n_trials) 
            
            #Update staircase randomness for next trial (https://questplus.readthedocs.io/en/latest/qp.html)
            SC.stim_selection_options["n"] += 7

            #Performance for current block
            num_correct = sum(acc) 
            tot_trials = len(acc)
            per_correct = sum(acc)/len(acc)
            mean_rt = np.nanmean(rt)
            #Initiate vectors new vectors 
            acc = [0] * n_trials
            rt = [0] * n_trials
                
            # offer a break
            points_text, speed_text, break_text, space, feedback_text = break_text_function(num_correct, tot_trials, mean_rt, per_correct, des_per_cor, des_mean_rt, blockN, n_blocks)
            points_text.draw(); speed_text.draw(); feedback_text.draw(); break_text.draw(); win.flip()
            core.wait(break_wait); points_text.draw(); speed_text.draw(); feedback_text.draw(); break_text.draw(); space.draw(); win.flip(); event.waitKeys(keyList=['space'])

            if blockN == scale_Nblock:
                ins.ExtraScale(win);core.wait(ins_wait); event.waitKeys(keyList=['lctrl']) # Press Control to continue the experiment (add in protocol)
            elif blockN == scale_Nblock + 1:
                ins.Main7(win);core.wait(ins_wait); event.waitKeys(keyList=['lctrl']) # Press Control to continue the experiment (add in protocol)
                
################################################################################################################################ 
# End of experiment
################################################################################################################################ 
break_text = vis.TextStim(win, text = "This is the end of the experiment. \n Thank you very much for your participation!")
space = vis.TextStim(win, text='Press space to close the experiment', pos=(0, -50), height=20)
break_text.draw(); win.flip()
core.wait(break_wait); break_text.draw(); space.draw(); win.flip(); event.waitKeys(keyList=['space'])
        
# Save data in a csv file -----------------------------------------------------
thisExp.saveAsWideText(file_name + '.csv', delim=',') 
  
edf_filename = f"Data/RDM_reportz_eyetrack_{sub}.EDF"
## save eyelink EDF from tracker to local Data folder
#tracker.setOfflineMode()
#tracker.closeDataFile()
try:
    print("Receiving EDF from EyeLink...")
    #tracker.receiveDataFile(edf_remote_name, edf_local_path)
    #print(f"EDF saved to {edf_local_path}")
except RuntimeError as e:
    print("Error transferring EDF:", e)

#tracker.close()

# End of the experiment -------------------------------------------------------
win.close()
core.quit() 


file_name = "Data/RDM_reportz_sub%d" % sub
thisExp = data.ExperimentHandler(dataFileName=file_name, extraInfo=info)  # saving extra info along with the main experimental data

