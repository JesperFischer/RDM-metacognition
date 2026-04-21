from psychopy import visual as vis
from psychopy import core, event
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

win = vis.Window(size=[1920,1080], units = 'pix', color='grey', allowGUI=False, fullscr=True)

def Intro(win):
    welcome_text = vis.TextStim(win, text = "Welcome to this experiment!", pos=(0,400), height = 55, color = 'white', wrapWidth=2000 )
    main_text = vis.TextStim(
        win,
        text=(
            "Today, you will be looking at moving dots and deciding whether you think\n"
            "they are moving mostly left or mostly right. We will practice this now.\n" 
            "Place your left ring finger on the “1” key, your left middle finger on the \n"
            "“2” key, and your left index finger on the  “3” key. Now place your right index \n"
            "finger on the “6” key, your right middle finger on the “7” key, and your \n"
            "right ring finger on the “8” key. If you think the dots are moving left, \n"
            "press the “3”, and if you think they are moving right, press the “6”\n"
            "You will get feedback about your decision.\n"
        ),
        pos=(0, 50) ,
        height= 50,
        color='white',
        wrapWidth=1700,
        alignText='left'
    )

    img = vis.ImageStim(
        win,
        image="keyboard_RDM_2.png",   
        pos=(400, -270),              
        size=(700, 300)                
    )
    
    hand_L = vis.ImageStim(
        win,
        image="hand_left_2.png",   
        pos=(0, -270),              
        size=(200,200)             
    )

    hand_R = vis.ImageStim(
        win,
        image="hand_right_2.png",   
        pos=(800, -270),              
        size=(200,200)             
    )
    
    space_text = vis.TextStim(
        win,
        text= "Press SPACE to continue.",
        pos=(0, -350) ,
        height= 35,
        color='white',
        wrapWidth=1400,
        alignText='left'
    )
 
    welcome_text.draw()
    main_text.draw()
    img.draw()
    
    hand_L.draw()
    hand_R.draw()
    space_text.draw() 
    win.flip()



def Tutorial(win):

    great = vis.TextStim(win, text = "Great Job!", pos = (0,350), height= 60,
        color='white',
        wrapWidth=1700,
        alignText='center')
    
    main_text = vis.TextStim(
        win,
        text=(
            "We’re now going to change the task a little. \n"
            "From now on, you will see the dots for a fixed amount of time.\n" 
            "When the dots disappear, that is your cue to respond.\n"
            "You must respond IMMEDIATELY after the dots disappear. \n"
            "If you wait too long, the trial will be missed. \n"
            "You cannot respond while the dots are still on the screen.\n"

            "Good Luck!"
        ),
        pos=(0, 50) ,
        height= 50,
        color='white',
        wrapWidth=1700,
        alignText='center'
    )
    space_text = vis.TextStim(
        win,
        text= "Press SPACE to continue.",
        pos=(0, -350) ,
        height= 35,
        color='white',
        wrapWidth=1400,
        alignText='center'
    )

    great.draw()
    main_text.draw()
    space_text.draw()
    win.flip()

def Main1(win):
    text_1 = vis.TextStim(win, text = "Good job! You seem to get the hang of the dots.", pos=(0,150), height = 40, color = 'white', wrapWidth=2000 )
    text_2 = vis.TextStim(win, text = "We will now introduce something new to the task.", pos=(0,0), height = 40, color = 'white', wrapWidth=2000 )

    space_text = vis.TextStim(
        win,
        text= "Call the experimenter.",
        pos=(0, -250) ,
        height= 30,
        color='white',
        wrapWidth=1400,
        alignText='center'
    )

    text_1.draw()
    text_2.draw()
    space_text.draw()
    win.flip()

def Main2(win):
    text_3 = vis.TextStim( win=win,
        text="In the next blocks, we will ask you to report the certainty around the direction of the dots. After your choice, you will see a scale like this:",
        pos=(0, 300),color='white',height= 40,wrapWidth=1200, alignText='center'
    )

    #confidence_question = vis.TextStim(win, text = "How confident were you in your decision?", pos=(0,100), height = 30, wrapWidth=1200)

    space_text = vis.TextStim(
        win,
        text= "Press SPACE to continue.",
        pos=(0, -325),
        height = 30,
        color='white',
        wrapWidth=1400,
        alignText='center'
    )

    slider = vis.Slider(win, ticks=[1, 2, 3, 4, 5, 6],   # 6 discrete positions
    labels=["Sure\n" "Left", "Left", "Slightly\n" "Left", "Slightly \n" "Right", "Right", "Sure\n" "Right"],
    granularity=1, style='rating',size = (1000,70), pos=(0, 100))


    #slider.marker.color = "white"
    #slider.markerPos = 0.5 
    #slider.marker.size = 30                       
    #slider_label_wrong = vis.TextStim(win, text= "definitely wrong", pos=(-300, 30), color = "white", height=23) 
    #slider_label_right = vis.TextStim(win, text= "definitely right", pos=(300, 30), color = "white", height=23)

    text_4 = vis.TextStim( win=win,
        text="We would like you to report your certainty in the direction. ",
        pos=(0, -175),color='white',height= 40,wrapWidth=1200, alignText='center'
    )

    text_3.draw()
    text_4.draw()
    space_text.draw()
    #confidence_question.draw()
    slider.draw()
    #slider_label_right.draw()
    #slider_label_wrong.draw()
    win.flip()

def Main3(win):
    text_5 = vis.TextStim( win=win,
        text="Each finger represents a specific direction and confidence level (as shown in the picture). " \
        "To indicate your answer, press the finger that best reflects your belief about the direction and your confidence in that decision.", 
       
        pos=(0, 400),color='white',height= 40,wrapWidth=1400, alignText='left'
    )

    slider = vis.Slider(win, ticks=[1, 2, 3, 4, 5, 6],   # 6 discrete positions
    labels=["Sure\n" "Left", "Left", "Slightly\n" "Left", "Slightly \n" "Right", "Right", "Sure\n" "Right"],
    granularity=1, style='rating',size = (1000,70), pos=(0, 200))

    img = vis.ImageStim(
        win,
        image="keyboard_RDM_2.png",   
        pos=(0, -150),              
        size=(700, 300)                
    )

    space_text = vis.TextStim(
        win,
        text= "Press SPACE to continue.",
        pos=(0, -350) ,
        height= 30,
        color='white',
        wrapWidth=1400,
        alignText='center'
    )

    hand_L = vis.ImageStim(
        win,
        image="hand_left_2.png",   
        pos=(-600, -150),              
        size=(300,300)             
    )

    hand_R = vis.ImageStim(
        win,
        image="hand_right_2.png",   
        pos=(600, -150),              
        size=(300,300)             
    )

    text_5.draw()
    slider.draw()
    img.draw()
    hand_L.draw()
    hand_R.draw()
    space_text.draw()
    win.flip()

Main2(win); event.waitKeys(keyList=['space']); win.close()