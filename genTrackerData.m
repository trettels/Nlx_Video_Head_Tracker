function [] = genTrackerData(indir, outfile, xmin, xmax, ymin, ymax, frameOffset)
	%% Function to estimate the head-direction of a rat connected to an HS-54 headstage, recorded on a Neuralynx system.
	%
	%	Input variables:
	%		indir: Directory containing the session data you are looking to analyze.  
	%			Must contain a video file named VT1.avi or VT0*.avi, a Neuralynx event
	%			file named Events.nev, and a neuralynx position file named VT1.nvt.
	%		outfile: The path and filename of the output file desired (will be saves as a *.mat file
	%		xmin: the left-most position x value of the recording arena
	%		xmax: the right-most position x value of the recording arena
	%		ymin: the bottom-most position y value of the recording arena
	%		ymax: the top-most position y value of the recording arena
	%		frameOffset: the time (in ms) of the first frame of this video, useful when multiple sessions are on one video.
	
	
	% Setup and initializations
	eval('cd ',indir);
	tic;
	global dir_delim;
	if ispc
		dir_delim='\';
	else
		dir_delim = '/';
	end
	printf('\n\n');
	% convert the MPEG-2 to an AVI
	fcheck = fopen('VT1.avi');
	if(fcheck < 0)
		trackvid = mpg2avi([]);
	else
		trackvid = VideoReader('VT1.avi');
		fclose(fcheck);
	end
	
	% Split the tracking epochs
	[nbegin, ~, t_begin, ~] = splitSessions('Events.nev');
	[timestamps]=readVideoTimestamps('VT1.nvt');
	
	
	% for each active trial...
	for aa=1:nbegin
		t_0=t_begin(aa,1);  %start time for epoch
		t_end = t_begin(aa,2); %end time for epoch
		f_0 = find(timestamps>=t_0, 1,'first')-frameOffset; %first frame in epoch
		f_end = find(timestamps<=t_end, 1, 'last')-frameOffset;  %last frame in epoch
		track = [];
		% for each frame within the trial...
		printf('\nTrial %i, number of frames: %i, t_0: %d, t_end: %d\n',aa, f_end-f_0, t_0, t_end);
		for bb= f_0:f_end
			if(mod(bb-f_0+1,1000)==0)
				fprintf('+');
			end
			if(mod(bb-f_0+1,10000)==0)
				fprintf('\n');
			end
			
			% Assign timepoint
			t_i=timestamps(bb+frameOffset);

			% Get the appropriate frame
			[rawframe] = read(trackvid, bb);
			scaledframe = max(0, double(rawframe)-repmat(sum(rawframe,3)/3 ,[1,1,3]) );
			correctedframe = removeReflections(scaledframe, xmin, xmax, ymin, ymax);
			redframe = correctedframe(:,:,1);
			greenframe = correctedframe(:,:,2);
			blueframe = correctedframe(:,:,3);
			
			% Determine x and y coordinates
			[x_i, y_i, rc, gc, ~]  = getHeadLocation(redframe, greenframe, blueframe);
			if(isempty(x_i) || isempty(y_i) )
				x_i=nan;
				y_i=nan;
				theta_i=nan;
				rc=[nan nan];
				gc=[nan nan];
			elseif(~isempty(rc) && ~isempty(gc))
				% Determine heading
				theta_i = getHeadAngle(rc, gc);
			else
				rc=[nan nan];
				gc=[nan nan];
				theta_i=nan;
			end
			%Assign the tracking info for the frame

			[ii, jj] = size(track);
			track=[track; t_i, x_i, y_i, theta_i, rc, gc]; % time, x, y, angle, red-center, green-center
		end
		% write tracking matrix to disk
		tk=toc;
		outdir = strcat(indir, dir_delim, 'begin',num2str(aa));
		exportTrackData(track, outdir, outfile);
	end
	printf('\nElapsed Time: %f seconds (%f minutes)',tk, tk/60);
	printf('\n\n');
	
	
	
%-----------------------------
% Helper Functions
%-----------------------------

%===Tracking info===%
%% remove reflections
function [frame] = removeReflections(frame, xmin, xmax, ymin, ymax)
	frame(1:ymin,:,:)=0; %top
	frame(ymax:end,:,:)=0; %bottom
	frame(:,1:xmin,:) = 0; %left
	frame(:,xmax:end,:) = 0; %right
	


%% Calculate the x and y coordinates
function [xx, yy, rc, gc, bc] = getHeadLocation(rf, gf, bf)
	
	[ri, rj]=find(rf>0);
	
	%keyboard
	if( ((length(ri) > 3)  && (length(rj)>3)) )
		nRed=collinearTest(ri,rj);
		if (nRed > 2 )   %ensure there are at least 3 non-collinear points
			rhull = [ri( convhull(ri, rj) ), rj( convhull(ri, rj) )];
			ry=real(centroid(rhull));
			rx=imag(centroid(rhull));
			rc = [rx, ry];
		else
			rc=[];
			rx=[];
			ry=[];
		end
	else
		rc=[];
		rx=[];
		ry=[];
	end
	
	[gi, gj]=find(gf>0);
	if( ((length(gi) > 3) && (length(gj)>3)) )
		nGreen=collinearTest(gi,gj);
		if(nGreen > 2) 
			ghull = [gi( convhull(gi, gj) ), gj( convhull(gi, gj) )];
			gy=real(centroid(ghull));
			gx=imag(centroid(ghull));
			gc=[gx, gy];
		else
			gc=[];
			gx=[];
			gy=[];
		end
	else
		gc=[];
		gx=[];
		gy=[];
	end
	
	[bi, bj]=find(bf>0);
	if( ((length(bi) > 3) && (length(bj)>3)) )
		nBlue=collinearTest(bi,bj);
		if(nBlue > 2)			
			bhull = [bi( convhull(bi, bj) ), bj( convhull(bi, bj) )];
			by=real(centroid(bhull));
			bx=imag(centroid(bhull));
			bc=[bx,by];
		else
			bc=[];
			bx=[];
			by=[];
		end
	else
		bc=[];
		bx=[];
		by=[];
	end
	
	
	if(isempty(gc) && isempty(rc))
		xx=[];
		yy=[];
	elseif(isempty(gc))
		xx=rx;
		yy=ry;
	elseif(isempty(rc))
		xx=gx;
		yy=gy;
	else
		xx=(rx+gx)/2;
		yy=(ry+gy)/2;
	end
	
	
function [nSlopes] = collinearTest(x,y)
	m=(x-x(1)) ./ (y-y(1));
	nSlopes=length(unique(m));
%% Calculate rat heading
function [ theta ] = getHeadAngle(rc, gc)	
	theta = atan2(rc(2)-gc(2),rc(1)-gc(1));
	

%===File Handling===%

%% use ffmpeg to convert mpeg-2 to avi/raw video
function [vidout] = mpg2avi(path)
	if(~isempty(path))
		eval('cd ', path);
	end
	cmd = ['ffmpeg -i VT1.mpg -vcodec rawvideo VT1.avi'];
	dos(cmd);
	vidout=VideoReader('VT1.avi');

%% load the event file and find the begin and sleep session timestamps
function [nbegin, nsleep, t_begin, t_sleep] = splitSessions(ev_file)
	FieldSelection = [1 0 0 0 1];  %Only get the timestamps and event strings
	ExtractHeader=0; %Don't get the header
	ExtractMode=1; %Extract all events
	ModeArray=[]; %not used
	[TimeStamps, EventStrings] = Nlx2MatEV( ev_file, FieldSelection, ExtractHeader, ExtractMode, ModeArray );
	nevents=length(EventStrings);
	
	nbegin=0;
	nsleep=0;
	t_begin=[];
	t_sleep=[];
	for aa=2:nevents-1
		switch lower(EventStrings{aa}(1:5))
			case('begin')
				nbegin = nbegin + 1;
				t_begin=[t_begin; TimeStamps(aa), TimeStamps(aa+1)]; %Double check that this is correct over [ts(aa), ts(aa+1)-1]
			case('sleep')
				nsleep = nsleep + 1;
				sleep=[t_sleep; TimeStamps(aa), TimeStamps(aa+1)];
			otherwise
				continue;
		end
	end
	
%% Load the full timestamps from the neuralynx tracking data
% Reads the video data from the video file (nvt) using the NeuraLynx dll
% for reading video data.
function posTime = readVideoTimestamps(posFile)
	% Want  timestamps and targets
	fieldSelect = [1,0,0,0,0,0];
	% Get header
	getHeader = 0;
	% Exctract every record
	extractMode = 1;

	% Get the data
	%[posData.t,posData.targets] = Nlx2MatVT_v4(posFile,fieldSelect,getHeader,extractMode);
	[posTime] = Nlx2MatVT(posFile,fieldSelect,getHeader,extractMode);

	
%% Write the new tracker info to file
function [] = exportTrackData(track, outdir, outfile)
	global dir_delim;
	filnam=strcat(outdir, dir_delim, outfile);
	save(filnam, 'track');