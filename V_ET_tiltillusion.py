#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy2 Experiment Builder (v1.85.4),
    on February 06, 2018, at 19:19
If you publish work using this script please cite the PsychoPy publications:
    Peirce, JW (2007) PsychoPy - Psychophysics software in Python.
        Journal of Neuroscience Methods, 162(1-2), 8-13.
    Peirce, JW (2009) Generating stimuli for neuroscience using PsychoPy.
        Frontiers in Neuroinformatics, 2:10. doi: 10.3389/neuro.11.010.2008
"""
from __future__ import absolute_import, division # must be in the first line
from psychopy import locale_setup, gui, monitors, visual, core, data, event, logging
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)
import numpy as np  # whole numpy lib is available, prepend 'np.'
import os  # handy system and path functions
import sys  # to get file system encoding
import psychopy.tools.monitorunittools as p2d

#To make the estimations of framerate more accurate
for frameN in range(200): 
    warmup=1
    
# store info about the experiment and experimental run
expName = 'V_tiltillusion'  # from the Builder filename that created this script
expInfo = {'participant':'', 'session':'001','run':'','ETused':[False,True]}
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName)
if dlg.OK == False:
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName


# %% setup for display
distanceMon = 57  # [58] in psychoph lab [99] in scanner
widthMon = 38  # [53] in psychoph lab [30] in scanner
PixW = 1280  # [1920.0] in psychopy lab [1920.0] in scanner
PixH = 1024  # [1080.0] in psychoph lab [1200.0] in scanner

moni = monitors.Monitor('testMonitor', width=widthMon, distance=distanceMon)
moni.setSizePix([PixW, PixH])  # [1920.0, 1080.0] in psychoph lab
expInfo['winsize'],expInfo['distanceMon'],expInfo['widthMon'] = str([PixW,PixH]),str(distanceMon),str(widthMon)

# %%Setup the Window
win = visual.Window(
    size=(PixW, PixH), screen=0,units = 'degFlat',
    allowGUI=False, allowStencil=False,
    monitor=moni, color=[0,0,0], colorSpace='rgb',
    blendMode='avg', useFBO=True)
# store frame rate of monitor if we can measure it
expInfo['frameRate'] = win.getActualFrameRate(nIdentical=60,nMaxFrames=100, nWarmUpFrames=10, threshold=1)
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess

#%% saving and logging

# get current path and save to variable _thisDir
# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__)).decode(sys.getfilesystemencoding())

# get parent path and move up one directory
strPathParentUp = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__),'..')),'data')
if not os.path.isdir(strPathParentUp):
    os.makedirs(strPathParentUp)
#move to parent_up_path
os.chdir(strPathParentUp)

# name and create specific subject folder
subjFolderName = strPathParentUp+os.path.sep+'%s_SubjData'%(expInfo['participant'])

if not os.path.isdir(subjFolderName):
    os.makedirs(subjFolderName)

#Name and create data folder for the experiment
dataFolderName = subjFolderName+os.path.sep+'%s'%(expInfo['expName'])
if not os.path.isdir(dataFolderName):
	os.makedirs(dataFolderName)
#Name and create specific folder for logging results
logFolderName = dataFolderName+os.path.sep+'Logging'
if not os.path.isdir(logFolderName):
	os.makedirs(logFolderName)
logFileName = logFolderName + os.path.sep + '%s_%s_Run%s_%s' %(expInfo['participant'], expInfo['expName'],
    expInfo['run'], expInfo['date'])

#Name and create specific folder for pickle output 
outFolderName = dataFolderName + os.path.sep + 'Output'
if not os.path.isdir(outFolderName):
	os.makedirs(outFolderName)

outFileName = outFolderName+os.path.sep+'%s_%s_Run%s_%s' %(expInfo['participant'], expInfo['expName'],
    expInfo['run'], expInfo['date'])


## Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
#filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath=None,
    savePickle=True, saveWideText=True,
    dataFileName=logFileName)

# save a log file for detail verbose info
logFile = logging.LogFile(logFileName+'.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file
# %%flag for 'escape' or other condition => quit the exp
endExpNow = False

# %%Start Code - component code to be run before the window creation
OriOfElems    = [90,105,105,105,105,105,105]
NumOfElements = np.array([1,6])
RadiusOfCircle = np.array([0.0,1])
PosOfEle = [[RadiusOfCircle[i]*np.cos(j*(2*np.pi/NumOfElements[i])),RadiusOfCircle[i]*np.sin(j*(2*np.pi/NumOfElements[i]))] for i in range(len(RadiusOfCircle)) for j in range(NumOfElements[i])]
PosOfEle = np.asarray(PosOfEle)
PosOffset = [4,0] # basic eccentricity


#def convert800600todeg(x,y): # convert eyelink coordinate sytem 800x600 (-400L:400R;-300U:300B) to visual degrees of the current monitor (via the pixel coordinate system of the current monitor)
#    #monsize = (800,600)
#    monsize = win.size
#    #monsize = win.monitor.getSizePix()
#    monsize = np.array(monsize)/2
#    x = p2d.pix2deg((x/400)*monsize[0],win.monitor)
#    y = p2d.pix2deg((y/300)*monsize[1],win.monitor)
#    return x,y
#
# generates conditions
Esize, Cpos, Espace, Cchg, CueFile, Sori = [1,1.5], [4,7], [1,1.5],[0,+1, +3, +7, -1, -3,-7],['<','>'],[75,105]
SF = 3
Tcond       = []
for n in range(len(Esize)):
    if Esize[n] == np.max(Esize):
        SF = 3
    else:
        SF = 3 * (Esize[1]/Esize[0])
    for m in range(len(Espace)):
        RadiusOfCircle = [0,Espace[m]*Esize[n]]
        p = [[RadiusOfCircle[g]*np.cos(h*(2*np.pi/NumOfElements[g])),RadiusOfCircle[g]*np.sin(h*(2*np.pi/NumOfElements[g]))] 
            for g in range(len(RadiusOfCircle)) for h in range(NumOfElements[g])]
        for i in range(len(CueFile)):
            if i == 0:                
                for j in range(len(Cchg)):
                    if Cchg[j] > 0:
                        corrans = 'up'
                    elif Cchg[j] < 0:
                        corrans = 'down'
                    else:
                        corrans = 'h'
                    for k in range(len(Sori)):
                        rori = [90]
                        lori = [90+Cchg[j]]
                        lori.extend([Sori[k] for s in range(6)])
                        rori.extend([Sori[k] for s in range(6)]) 
                        Tcond.append({'esize':Esize[n],'pos':p,'Cue':CueFile[i], 'cchg':Cchg[j],'lori':lori,'rori':rori,'corrans':corrans,'sf':SF})

            if i == 1:
                
                for j in range(len(Cchg)):
                    if Cchg[j] > 0:
                        corrans = 'right'
                    elif Cchg[j] < 0:
                        corrans = 'left'
                    else:
                        corrans = 'Vertical'
                    for k in range(len(Sori)):
                        lori = [90]
                        rori = [90+Cchg[j]]
                        rori.extend([Sori[k] for s in range(6)])
                        lori.extend([Sori[k] for s in range(6)])                
                        Tcond.append({'esize':Esize[n],'pos':p,'Cue':CueFile[i], 'cchg':Cchg[j],'lori':lori,'rori':rori,'corrans':corrans,'sf':SF})
                    
np.random.shuffle(Tcond)
# %% define functions
#def convert800600todeg(x,y): # convert eyelink coordinate sytem 800x600 (-400L:400R;-300U:300B) to visual degrees of the current monitor (via the pixel coordinate system of the current monitor)
#    #monsize = (800,600)
#    monsize = win.size
#    #monsize = win.monitor.getSizePix()
#    monsize = np.array(monsize)/2
#    x = p2d.pix2deg(x-(monsize[0]/2),win.monitor)
#    y = p2d.pix2deg(x-(monsize[0]/2),win.monitor)
#    return x,y
# %%
import eyeTracker

# %%Initialize components for Routine "BeginIns"
BeginInsClock = core.Clock()
text_expin = visual.TextStim(win=win, name='expins',
    text='Please focus on the center cross and cue.\n\nPress up key when you detect an anticlockwise center grating. \n\n Press down key when you detect an clockwise center grating. \n\n',
    pos=(0, 0), height=1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1,
    depth=0.0);

# %%Initialize components for Routine "trial"
trialClock = core.Clock()


fix = visual.TextStim(win=win, name='fix',
    text='default text',
    font='Arial',
    pos=(0, 0), height=1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1,
    depth=0.0);
cue = visual.TextStim(win=win, name='cue',
    text='default text',
    font=u'Arial',
    pos=(0, 0), height=1, wrapWidth=None, ori=0, 
    color=u'white', colorSpace='rgb', opacity=1,
    depth=-1.0);
Lgrating = visual.ElementArrayStim(win=win, 
    elementMask='circle', oris = OriOfElems,
    nElements=int(np.sum(NumOfElements)), 
    sizes=1, sfs=3, contrs=0.7, xys= PosOfEle - PosOffset, 
    colors=(1, 1, 1), colorSpace='rgb')
Rgrating = visual.ElementArrayStim(win=win,
    elementMask='circle', oris = OriOfElems,
    nElements=int(np.sum(NumOfElements)), 
    sizes=1, sfs=3, contrs=0.7, xys= PosOfEle + PosOffset, 
    colors=(1, 1, 1), colorSpace='rgb')
text = visual.TextStim(win=win, name='text',
    text=u'',
    font=u'Arial',
    pos=(0, 0), height=1, wrapWidth=None, ori=0, 
    color=u'white', colorSpace='rgb', opacity=1,
    depth=-4.0);
# %%Initialize components for Routine "rest"
restClock = core.Clock()
text_rest = visual.TextStim(win=win, name='text_rest',
    text=u'Have a Rest \n\nPress any key to start',
    font=u'Arial',
    pos=(0, 0), height=1, wrapWidth=None, ori=0, 
    color=u'white', colorSpace='rgb', opacity=1,
    depth=0.0);


# %%Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 
# %%Initiate the eyetracker
eyetracker = eyeTracker.Tracker_EyeLink(win=win, clock=globalClock,sj=expInfo['participant']) # enable when using eyelink 1000
#eyetracker = eyeTracker.Tracker_Dummy(win=win, clock=globalClock,sj=expInfo['participant']) # enable when using mouse simulation for eye_movements
#Sent the initial setup to the eyetracker
st = "eyexp_start_" + expInfo['participant']
eyetracker.sendMessage(st)    # the original function is just to print the content of st
core.wait(1)

## ------Prepare to start Routine "expin"-------
t = 0
BeginInsClock.reset()  # clock
frameN = -1
continueRoutine = True
# update component parameters for each repeat
key_resp_expins = event.BuilderKeyResponse()
# keep track of which components have finished
expinComponents = [text_expin, key_resp_expins]
for thisComponent in expinComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED

# -------Start Routine "expin"-------
eyetracker.sendMessage('Expinstruction')
while continueRoutine:
    # get current time
    t = BeginInsClock.getTime()
#    eyetracker.preTrial(trial=t,calibTrial=0,win=win) #eyelink drift correction
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *text_expin* updates
    if t >= 0.0 and text_expin.status == NOT_STARTED:
        # keep track of start time/frame for later
        text_expin.tStart = t
        text_expin.frameNStart = frameN  # exact frame index
        text_expin.setAutoDraw(True)
#        eyetracker.sendMessage('expinsF')
    #    frameRemains = 0.0 + 1.0- win.monitorFramePeriod * 0.75  # most of one frame period left
    #    if text_expin.status == STARTED and t >= frameRemains:
    if text_expin.status == STARTED :
        text_expin.setAutoDraw(True)
    
    # *key_resp_expins* updates
    if t >= 0.0 and key_resp_expins.status == NOT_STARTED:
        # keep track of start time/frame for later
        key_resp_expins.tStart = t
        key_resp_expins.frameNStart = frameN  # exact frame index
        key_resp_expins.status = STARTED
#        eyetracker.sendMessage('expins K    ')
        # keyboard checking is just starting
        win.callOnFlip(key_resp_expins.clock.reset)  # t=0 on next screen flip
        event.clearEvents(eventType='keyboard')
    if key_resp_expins.status == STARTED:
        theseKeys = event.getKeys()
        
        # check for quit:
        if "escape" in theseKeys:
            endExpNow = True

        
        if len(theseKeys) > 0:  # at least one key was pressed
            key_resp_expins.keys = theseKeys[-1]  # just the last key pressed
            key_resp_expins.rt = key_resp_expins.clock.getTime()
            # a response ends the routine
            continueRoutine = False
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in expinComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # check for quit (the Esc key)
    if endExpNow or event.getKeys(keyList=["escape"]):
        eyetracker.sendMessage('endinterrupted')
        core.quit()
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "expin"-------
for thisComponent in expinComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# check responses
if key_resp_expins.keys in ['', [], None]:  # No response was made
    key_resp_expins.keys=None
thisExp.addData('key_resp_expins.keys',key_resp_expins.keys)
if key_resp_expins.keys != None:  # we had a response
    thisExp.addData('key_resp_expins.rt', key_resp_expins.rt)
thisExp.nextEntry()
# the Routine "expin" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# set up handler to look after randomisation of conditions etc
block = data.TrialHandler(nReps=10, method='random',
                          extraInfo=expInfo, originPath=-1,
                          trialList=[None],
                          seed=None, name='block')
thisExp.addLoop(block)  # add the loop to the experiment
thisBlock = block.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisBlock.rgb)
if thisBlock != None:
    for paramName in thisBlock.keys():
        exec(paramName + '= thisBlock.' + paramName)

for thisBlock in block:
    currentLoop = block
    # abbreviate parameter names if possible (e.g. rgb = thisBlock.rgb)
    if thisBlock != None:
        for paramName in thisBlock.keys():
            exec(paramName + '= thisBlock.' + paramName)

        # set up handler to look after randomisation of conditions etc
    trials = data.TrialHandler(nReps=2, method='random',
                                   extraInfo=expInfo, originPath=-1,
                                   trialList=Tcond,
                                   seed=None, name='trials')
    thisExp.addLoop(trials)  # add the loop to the experiment
    thisTrial = trials.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
    if thisTrial != None:  
        for paramName in thisTrial.keys():
            exec(paramName + '= thisTrial.' + paramName)
    
        for thisTrial in trials:
            currentLoop = trials
            eyetracker.preTrial(trial=currentLoop.thisTrialN,calibTrial=0,win=win) #eyelink drift correction
            # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
            if thisTrial != None:
                for paramName in thisTrial.keys():
                    exec(paramName + '= thisTrial.' + paramName)
            
            # ------Prepare to start Routine "trial"-------
            t = 0
            fixroi,eye_pos = [[]for i in range(2)]
            trialClock.reset()  # clock
            frameN = -1
            continueRoutine = True
            routineTimer.add(10)
            # update component parameters for each repeat
            fix.setText('+')
            cue.setText(Cue)
            Lgrating.setXYs(np.asarray(pos)- PosOffset)
            Lgrating.setOris(lori)
            Lgrating.setSizes(esize)
            Lgrating.setSfs(sf)
            Rgrating.setXYs(np.asarray(pos)+ PosOffset)
            Rgrating.setOris(rori)
            Rgrating.setSizes(esize)
            Rgrating.setSfs(sf)
            key_resp_2 = event.BuilderKeyResponse()
            # keep track of which components have finished
            trialComponents = [fix, cue, Lgrating, Rgrating, text, key_resp_2]
            for thisComponent in trialComponents:
                if hasattr(thisComponent, 'status'):
                    thisComponent.status = NOT_STARTED
    
            # -------Start Routine "trial"-------
            #        while continueRoutine and routineTimer.getTime() > 0:
            currenttrial = str(currentLoop.thisTrialN) + '_' + str(currentLoop.thisRepN)
            eyetracker.sendMessage(currenttrial)
            while continueRoutine :
                # get current time
                t = trialClock.getTime()
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                # *fix* updates
                if t >= 0.0 and fix.status == NOT_STARTED:
                    # keep track of start time/frame for later
                    fix.tStart = t
                    fix.frameNStart = frameN  # exact frame index
                    fix.setAutoDraw(True)
                frameRemains = 0.0 + 0.3- win.monitorFramePeriod * 0.75  # most of one frame period left
                if fix.status == STARTED and t >= frameRemains:
                    fix.setAutoDraw(False)
            
                # *cue* updates
                if t >= 0.3 and cue.status == NOT_STARTED:
                # keep track of start time/frame for later
                    cue.tStart = t
                    cue.frameNStart = frameN  # exact frame index
                    cue.setAutoDraw(True)
                frameRemains = 0.3 + 0.5- win.monitorFramePeriod * 0.75  # most of one frame period left
                if cue.status == STARTED and t >= frameRemains:
                    cue.setAutoDraw(False)
    
                # *Lgrating* updates
                if t >= 0.6 and Lgrating.status == NOT_STARTED:
                    # keep track of start time/frame for later
                    Lgrating.tStart = t
                    Lgrating.frameNStart = frameN  # exact frame index
                    Lgrating.setAutoDraw(True)
                frameRemains = 0.6 + 0.2- win.monitorFramePeriod * 0.75  # most of one frame period left
                if Lgrating.status == STARTED and t >= frameRemains:
                    Lgrating.setAutoDraw(False)
        
                # *Rgrating* updates
                if t >= 0.6 and Rgrating.status == NOT_STARTED:
                # keep track of start time/frame for later
                    Rgrating.tStart = t
                    Rgrating.frameNStart = frameN  # exact frame index
                    Rgrating.setAutoDraw(True)
                frameRemains = 0.6 + 0.2- win.monitorFramePeriod * 0.75  # most of one frame period left
                if Rgrating.status == STARTED and t >= frameRemains:
                    Rgrating.setAutoDraw(False)
            
                # *text* updates
                if t >= 0.8 and text.status == NOT_STARTED:
                    # keep track of start time/frame for later
                    text.tStart = t
                    text.frameNStart = frameN  # exact frame index
                    text.setAutoDraw(True)
                frameRemains = 0.8 + 0.8- win.monitorFramePeriod * 0.75  # most of one frame period left
                if text.status == STARTED and t >= frameRemains:
                    text.setAutoDraw(False)
    
                # *key_resp_2* updates
                if t >= .9 and key_resp_2.status == NOT_STARTED:
                    # keep track of start time/frame for later
                    key_resp_2.tStart = t
                    key_resp_2.frameNStart = frameN  # exact frame index
                    key_resp_2.status = STARTED
    #                eyetracker.sendMessage('expK')
                    # keyboard checki ng is just starting
                    win.callOnFlip(key_resp_2.clock.reset)  # t=0 on next screen flip
                    event.clearEvents(eventType='keyboard')
                #            frameRemains = 0.8 + 0.8- win.monitorFramePeriod * 0.75  # most of one frame period left
                #            if key_resp_2.status == STARTED and t >= frameRemains:
                #                key_resp_2.status = STOPPED
                #            if key_resp_2.status == STARTED:
                #                theseKeys = event.getKeys(keyList=['right', 'left'])
                            #Mark different period
                #record eyemovement
                eye_pos_int = eyetracker.getSample()
                xg,yg   = p2d.pix2deg(eye_pos_int[0],win.monitor),p2d.pix2deg(eye_pos_int[1],win.monitor)
                eye_pos.append([xg,yg])
                eye_shift = np.sqrt(xg**2 + yg**2)
                if eye_shift < 5  : # if eyemovements is larger than 0.3 degree, go to next trial
                    fixroi.append(1)
                else:
                    
                    fixroi.append(0)
                    # a response ends the routine
                    #continueRoutine = False
                if key_resp_2.status == STARTED:
                    theseKeys = event.getKeys(keyList=['up', 'down'])
                    # check for quit:
                    if "escape" in theseKeys:
                        endExpNow = True
                    if len(theseKeys) > 0:  # at least one key was pressed
                        key_resp_2.keys.extend(theseKeys)  # storing all keys
                        key_resp_2.rt.append(key_resp_2.clock.getTime())
                        # was this 'correct'?
                        if (key_resp_2.keys == [str(corrans)]) or (key_resp_2.keys == corrans):
                            key_resp_2.corr = 1
                        else:
                            key_resp_2.corr = 0
                        # a response ends the routine
                        continueRoutine = False
                # check if all components have finished
                if not continueRoutine:  # a component has requested a forced-end of Routine
                    break
                continueRoutine = False  # will revert to True if at least one component still running
                for thisComponent in trialComponents:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
        
                # check for quit (the Esc key)
                if endExpNow or event.getKeys(keyList=["escape"]):
    #                eyetracker.sendMessage('endinterrupted')
                    core.quit()
                
                # refresh the screen
                if continueRoutine:
    #                core.wait(0.5)# don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            core.wait(np.random.random())
    
            # -------Ending Routine "trial"-------
            for thisComponent in trialComponents:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            # check responses
            if key_resp_2.keys in ['', [], None]:  # No response was made
                key_resp_2.keys=None
                if str(corrans).lower() == 'none':
                    key_resp_2.corr = 1  # correct non-response
                else:
                    key_resp_2.corr = 0  # failed to respond (incorrectly)
    
    
        # store data for trials (TrialHandler)
            trials.addData('key_resp_2.rt',key_resp_2.rt)
            trials.addData('key_resp_2.keys',key_resp_2.keys)
            trials.addData('key_resp_2.corr', key_resp_2.corr)
            trials.addData('eyepos',eye_pos)
            trials.addData('fixroi',fixroi)
            if key_resp_2.keys != None:  # we had a response
                trials.addData('key_resp_2.rt', key_resp_2.rt)
            thisExp.nextEntry()
    
        # completed 5 repeats of 'trials'
    

    # ------Prepare to start Routine "rest"-------
    t = 0
    restClock.reset()  # clock
    frameN = -1
    continueRoutine = True
    # update component parameters for each repeat
    key_resp_rest = event.BuilderKeyResponse()
    # keep track of which components have finished
    restComponents = [text_rest, key_resp_rest]
    for thisComponent in restComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED

    # -------Start Routine "rest"-------
    eyetracker.sendMessage('restFix')
    while continueRoutine:
        # get current time
        t = restClock.getTime()
#        eyetracker.preTrial(trial=t,calibTrial=0,win=win) #eyelink drift correction
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *text_rest* updates
        if t >= 0.0 and text_rest.status == NOT_STARTED:
            # keep track of start time/frame for later
            text_rest.tStart = t
            text_rest.frameNStart = frameN  # exact frame index
            text_rest.setAutoDraw(True)
#            eyetracker.sendMessage('restFix')

        # *key_resp_rest* updates
        if t >= 0.0 and key_resp_rest.status == NOT_STARTED:
        # keep track of start time/frame for later
            key_resp_rest.tStart = t
            key_resp_rest.frameNStart = frameN  # exact frame index
            key_resp_rest.status = STARTED
#            eyetracker.sendMessage('restK')
            # keyboard checking is just starting
            win.callOnFlip(key_resp_rest.clock.reset)  # t=0 on next screen flip
            event.clearEvents(eventType='keyboard')
        if key_resp_rest.status == STARTED:
            theseKeys = event.getKeys()
            
            # check for quit:
            if "escape" in theseKeys:
#                eyetracker.sendMessage('endinterrupted')
                endExpNow = True
            if len(theseKeys) > 0:  # at least one key was pressed
                key_resp_rest.keys = theseKeys[-1]  # just the last key pressed
                key_resp_rest.rt = key_resp_rest.clock.getTime()
                # a response ends the routine
                continueRoutine = False

        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in restComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished

        # check for quit (the Esc key)
        if endExpNow or event.getKeys(keyList=["escape"]):
  #            eyetracker.sendMessage('endinterrupted')
            core.quit()
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()

    # -------Ending Routine "rest"-------
    for thisComponent in restComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # check responses
    if key_resp_rest.keys in ['', [], None]:  # No response was made
        key_resp_rest.keys=None
    block.addData('key_resp_rest.keys',key_resp_rest.keys)
    if key_resp_rest.keys != None:  # we had a response
        block.addData('key_resp_rest.rt', key_resp_rest.rt)
    # the Routine "rest" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    thisExp.nextEntry()

# completed 4 repeats of 'block'

# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(logFileName+'.csv')
thisExp.saveAsPickle(logFileName)
logging.flush()
# make sure everything is closed down
thisExp.abort()  # or data files will save again on exit
win.close()
eyetracker.closeConnection()
core.quit()