% This program tests heading discrimination for flow fields made of
% plaids. The envelopes are stationary, the patterns move. We're using
% up/down staircases to find observers' subjective forwards. Start
% values are randomly selected from the range of legal values.
%
%
% DMG
% Last Edited: 04/07/2015

clear all;close all;clc;ListenChar(2);HideCursor;                           % Start with a blank slate
try
    %-----Subject Settings------------------
    initials                        = 'test';
    n_trials                        = 25;                                   % Number of trials per staircase
    n_staircases                    = 2;                                    % Per condition
    
    stimulus_duration               = 3;                                   % secs
    
    translation_speed               = 2;                                    % m/s, apparently 1.9 m/s is a brisk walking speed
    rotation                        = [-1 0 1];                             % deg/s
    
    %-----Experiment Settings, Don't Change----
    %-----Movement Directions
    trans_axes                      = [3 1];                                % which way are we driving? 1=X,2=Y,3=Z first number is main direction, second is judged direction
    rot_axis                        = [0 1 0];                              % X (pitch),Y (yaw),Z (roll)
    
    %-----Array Settings
    plane_dist                      = [12.5 25];                            % m
    dot_density                     = [.16 .16];                            % dots/deg^2 on screen
    
    %-----Staircase Settings
    step                            = [4 2 1];
    limits                          = [-20 20];
    
    %-----Element Settings
    gabor_diam                      = 30;                                   % Arcmin
    gabor_sf                        = 2;                                    % c/deg
    gabor_contrast                  = 100;                                  % percent, for the pattern; each element will be half this
    spatial_envelope                = 2;                                    % 0 = disk, 1 = Gabor, 2 = raised cosine
    
    view_window                     = [36 27];                              % X,Y centered around fixation, in dva
    exclude                         = [-20 -3 20 3];                        % Don't place elements where the FOE will be
    
    predisplay                      = 0;                                    % show stationary world before it moves (sec)
    
    H_ecc_fix                       = 0;                                    % Horizontal fixation ecc (degs, neg is left)
    V_ecc_fix                       = 0;                                    % Vertical fixation ecc (degs, neg is up)
    
    experiment_id                   = 'pattern_heading';                         % Used in group filename
    fast                            = 1;                                    % Automatically trigger trials
    ITI                             = 2;                                    % Intertrial Inverval, Seconds
    background                      = 127;                                  % Grayscale Units
    fixate                          = 1;                                    % Present fixation spot during motion
    data_path                       = '/Users/heegerlab/Documents/Heading/Data/';   % Folder for saving data files
    %-----Rig Settings----------------------
    view_dist                       = .57;                                  % m
    scale_factor                    = 1.78;                                 % Arcmin/pixel
    frame_rate                      = 120;                                   % Screen frame rate (hz)
    linearize                       = 0;                                    % Use calibrated LUT (do this when available)
    gamma_fname                     = 'TMSgamma_042114';                    % Calibration table
    
    %-----Housekeeping----------------------
    % Scale things based on viewing distance, and convert other stuff to
    % the units PsychToolbox wants...
    nPlanes                         = length(plane_dist);
    tme                             = clock;
    logname = strcat(data_path,initials,'_',experiment_id,'_trials');
    
    gabor_size                      = gabor_diam/scale_factor;
    stimulus_radius                 = round(gabor_size/2);
    H_ecc_fix                       = H_ecc_fix*60/scale_factor;
    V_ecc_fix                       = V_ecc_fix*60/scale_factor;
    mv_length                       = ceil(stimulus_duration*frame_rate);
    f                               = (gabor_sf*scale_factor/60)*2*pi;
    angle                           = 0;
    a                               = cos(angle)*f;
    b                               = sin(angle)*f;
    amplitude                       = background;
    
    %-----Spatial Envelope------------------
    [x,y]=meshgrid(-stimulus_radius:stimulus_radius,-stimulus_radius:stimulus_radius);
    bps = (stimulus_radius)*2+1;
    circle=((stimulus_radius)^2-(x.^2+y.^2));
    for i=1:bps; for j =1:bps; if circle(i,j) < 0; circle(i,j) = 0; else circle(i,j) = 1; end; end;
    end;
    if spatial_envelope == 1
        circle = (exp(-(((x)/(sqrt(2)*Gaussian_stdev/6)).^2)-((y/(sqrt(2)*Gaussian_stdev/2)).^2)).*circle);
    elseif spatial_envelope == 2
        R = (sqrt(x.^2 + y.^2) + eps).*circle;
        R = R/max(max(R));
        cos2D = (cos(R*pi)+1)/2;
        circle = (cos2D.*circle);
    end
    circle = circle*255/2;
    
    %-----Open Screens----------------------
    oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
    oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
    screens=Screen('Screens');
    screenNumber=max(screens);
    
    w=Screen('OpenWindow',screenNumber,background,[],[],2);
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    screen_rect = Screen('Rect',w);
    %     if size(screens,2) == 2;
    %         w2=Screen('OpenWindow',0,0,[],[],2);
    %         Screen('FillRect',w2, 0); Screen('Flip', w2);
    %     end
    if linearize
        fid = fopen(gamma_fname,'r');
        screen_clut = fread(fid,[256 3],'float64');
        fclose(fid);
        screen_clut =   screen_clut -1;
        screen_clut = screen_clut/255;
        screen('LoadNormalizedGammaTable',screenNumber,screen_clut);
    end
    Screen('FillRect',w, background);
    Screen('Flip', w);
    Screen('FillRect',w, background);
    Screen('TextSize',w,20);Screen('TextFont',w,'Charcoal');
    
    %-----Screen Landmarks------------
    sr_hor = round(screen_rect(3)/2); % Middle of the screen, horizontally, in pixels
    sr_ver = round(screen_rect(4)/2); % Middle of the screen, vertically, in pixels
    fix_hor = sr_hor+H_ecc_fix;     % Horizontal location of fixation cross, in pixels
    fix_ver = sr_ver+V_ecc_fix;     % Vertical location of fixation cross, in pixels
    movie_rect= [0,0,bps,bps];
    
    %-----Make the grating in all phases-----
    gabor = zeros(1,360);
    for i=1:360
        grating = round(((sin(a*x+b*y+i*pi/180)*amplitude)+background));
        gabor(i) = Screen('MakeTexture',w,cat(3,grating,circle));
    end
    
    
    %-----Set Up Conditions-----------
    n_conditions = 0;
    for j=1:length(rotation)
        for k=1:n_staircases
            n_conditions = n_conditions+1;
            
            cond(n_conditions).rot_index = j;
            
            cond(n_conditions).rotate = [0 0 0];
            cond(n_conditions).rotate = rot_axis*rotation(j)*pi/180;
            
            cond(n_conditions).contrast = gabor_contrast/100;
            cond(n_conditions).count = 0;
            
            cond(n_conditions).last_resp = 0;
            cond(n_conditions).n_flips = 0;
            
            cond(n_conditions).heading = limits(1)+ceil(rand*(limits(2)-limits(1)));
            cond(n_conditions).step = step(1);
            
            cond(n_conditions).resp_history = zeros(1,n_trials);
            cond(n_conditions).heading_history = zeros(1,n_trials);
        end
    end
    results = zeros(n_conditions,n_trials);
    
    %-----Randomize Trials------------
    total_trials = n_trials*n_conditions;
    perm = randperm(total_trials);
    perm = mod(perm,n_conditions)+1;
    duration_check = zeros(size(perm));
    
    Screen('DrawText',w,'Heading discrimination: Pattern',100,100,250);
    Screen('DrawText',w,'Use Left/Right Arrows to respond',100,130,250);
    Screen('DrawText',w,'press SPACE BAR to start',100,160,250);
    Screen('Flip',w);
    
    FlushEvents('keyDown');
    validKey = 0;
    while ~validKey
        [secs, keyCode, deltaSecs] = KbWait(-1);
        if keyCode(KbName('space'))
            validKey = 1;
        end
    end
    
    Screen('FillRect', w, background);
    Screen('Flip', w);
    tic;
    %-----Main experimental loop-----------------
    for trial=1:total_trials
        aa = GetSecs;
        % Draw the white fixation cross
        Screen('FillRect',w, background);
        Screen('DrawDots', w, [fix_hor; fix_ver], 4, 250, [], 2);
        Screen('Flip',w);
        
        % Set up the world
        translate = [0 0 0];
        translate(trans_axes(1)) = translation_speed*cosd(cond(perm(trial)).heading);
        translate(trans_axes(2)) = translation_speed*sind(cond(perm(trial)).heading);
        
        angle_all = [];
        speed_all = [];
        pos_all = [];
        for i=1:nPlanes
            [dots_m, dots_deg] = make_dot_plane(dot_density(i), plane_dist(i), view_window, exclude);
            velocity_field = calculate_plane_flow(dots_m, dots_deg, translate, cond(perm(trial)).rotate, view_dist);
            nGabors(i) = size(dots_deg,2);
            posA = CenterRectOnPoint(movie_rect,dots_deg(1,:)'*60/scale_factor+sr_hor,dots_deg(2,:)'*60/scale_factor+sr_ver)'; %big Y numbers are lower on the screen in pixels so we have to flip
            pos = [posA posA];
            
            speed_ver = velocity_field(2,:);
            speed_hor = velocity_field(1,:);
            
            angle_ver = 270*(speed_ver>=0)+90*(speed_ver<0); % in deg
            angle_hor = 180*(speed_hor>0); % in deg
            speed_ver = abs(velocity_field(2,:));
            speed_hor = abs(velocity_field(1,:));
            
            speed = zeros(1,length(speed_ver)*2);
            angle = speed;
            speed = [speed_hor speed_ver];
            angle = [angle_hor angle_ver];
            angle_all = [angle_all angle];
            speed_all = [speed_all speed];
            pos_all = [pos_all pos];
        end
        angle = angle_all;
        speed = speed_all;
        pos = pos_all;
        phase_step = round(speed*gabor_sf*360/frame_rate)';
        
        cond(perm(trial)).count = cond(perm(trial)).count+1;
        phases = repmat(ceil(rand(size(speed,2),1)*360),1,mv_length)+repmat(phase_step,1,mv_length).*repmat(0:(mv_length-1),size(speed,2),1);
        phases = rem(phases,360);
        phases = phases + (phases==0)*360;
        
        % % Set Correct Answer
        % if translate(trans_axes(2))<0
        %     correct = 'LeftArrow';
        %     incorrect = 'RightArrow';
        % else
        %     correct = 'RightArrow';
        %     incorrect = 'LeftArrow';
        % end
        
        %Finish the ITI
        WaitSecs(ITI-(GetSecs-aa));
        
        % Draw the red fixation cross
        Screen('DrawDots', w, [fix_hor; fix_ver], 4, [250 0 0], [], 2);
        Screen('Flip',w);
        
        FlushEvents('keyDown');
        priorityLevel=MaxPriority(w);
        Priority(priorityLevel);
        if predisplay
            Screen('DrawTextures', w, gabor(phases(:,1)),movie_rect,pos,angle);
            if fixate
                Screen('DrawDots', w, [fix_hor; fix_ver], 4, [250 0 0], [], 2);
            end
            Screen('Flip',w);
            WaitSecs(predisplay);
        end
        
        % Play the movie
        StimulusOnsetTime = zeros(1,mv_length);
        aa = GetSecs;
        tag = clock;
        for frame = 1:mv_length
            Screen('DrawTextures', w, gabor(phases(:,frame)),movie_rect,pos,angle);
            if fixate
                Screen('DrawDots', w, [fix_hor; fix_ver], 4, [250 0 0], [], 2);
            end
            [VBLTimestamp, StimulusOnsetTime(frame), FlipTimestamp, Missed, Beampos] = Screen('Flip',w);
        end
        duration_check(trial) = GetSecs-aa;
        Screen('FillRect',w, background);
        if fixate
            Screen('DrawDots', w, [fix_hor; fix_ver], 4, [250 0 0], [], 2);
        end
        Screen('Flip',w);
        % Get the response
        validKey = 0;
        while ~validKey
            [secs, keyCode, deltaSecs] = KbWait(-1);
            if keyCode(KbName('ESCAPE'))
                Screen('CloseAll');
                Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
                Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
                ListenChar(1);
                ShowCursor;
                break
            elseif keyCode(KbName('LeftArrow'))
                validKey = 1;
                resp = 0;
            elseif keyCode(KbName('RightArrow'))
                validKey = 1;
                resp = 1;
                results(perm(trial),cond(perm(trial)).count) = 1;
                % elseif keyCode(KbName(incorrect))
                %     validKey = 1;
                % elseif keyCode(KbName(correct))
                %     validKey = 1;
                %     results(perm(trial),cond(perm(trial)).count) = 1;
            else
                Beeper('low');
            end
        end
        Priority(0);
        % Log the trial @@@
        fid = fopen(logname,'a');
        fprintf(fid,'%s, %s, %d, %d, %d, %d, %d, %f, %d, %f, %f, %f, %d, %f, %f\n', initials, experiment_id, tag(1), tag(2), tag(3), tag(4), tag(5), tag(6), trial, translation_speed, rotation(cond(perm(trial)).rot_index), cond(perm(trial)).heading, resp, stimulus_duration, duration_check(trial));
        fclose(fid);
        % Tell the staircase what happened
        cond(perm(trial)).resp_history(cond(perm(trial)).count) = resp;
        cond(perm(trial)).heading_history(cond(perm(trial)).count) = cond(perm(trial)).heading;
        % Was this a flip? If so, decrease the step size
        if resp ~= cond(perm(trial)).last_resp && cond(perm(trial)).count > 1
            cond(perm(trial)).n_flips = cond(perm(trial)).n_flips + 1;
            if cond(perm(trial)).n_flips < length(step)
                cond(perm(trial)).step = step(cond(perm(trial)).n_flips+1);
            else
                cond(perm(trial)).step = step(length(step));
            end
        end
        % Set the new heading
        if resp
            cond(perm(trial)).heading = cond(perm(trial)).heading-cond(perm(trial)).step;
            if cond(perm(trial)).heading < limits(1)
                cond(perm(trial)).heading = limits(1);
            end
        else
            cond(perm(trial)).heading = cond(perm(trial)).heading+cond(perm(trial)).step;
            if cond(perm(trial)).heading > limits(2)
                cond(perm(trial)).heading = limits(2);
            end
        end
        
        % Set the new last value
        cond(perm(trial)).last_resp = resp;
        
    end
    Screen('CloseAll');
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    ListenChar(1);
    clc;
    % Print out results here.
    for i=1:n_conditions
       if rotation(cond(i).rot_index) > 0
           plot(cond(i).heading_history,'bo-');
           hold on;
       elseif rotation(cond(i).rot_index) < 0
           plot(cond(i).heading_history,'ro-');
           hold on;
       else
           plot(cond(i).heading_history,'ko-');
           hold on;
       end
    end
    filename = strcat(data_path,initials,'_',experiment_id,'_',int2str(tme(1)),'_',int2str(tme(2)),'_',int2str(tme(3)),'_',int2str(tme(4)),'_',int2str(tme(5)));
    fprintf('\n');
    time = toc/60;
    fprintf('Elapsed time  (minutes)  =  %4.1f\n',time);
    clear x y R grating cos2D circle
    save(filename);
catch ME
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    ListenChar(1);
    ShowCursor;
    Screen('CloseAll');
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    Priority(0);
end %try..catch..
