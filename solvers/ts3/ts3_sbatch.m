function job = ts3_sbatch(job)
%QSUB Creates an MPI submission script for filename and submits it to the queue
%   The nodes argument is optional and is assumed to be 1 unless provided
%
%   The ppn argument is optional and is assumed to be 4 unless provided
%
%   The walltime argument is optional and is assumed to be 12:00:00
%   unless provided. If provided, it is given as a floating point number of
%   hours, eg. 2.5 for 2 hours, 30 minutes.
%
%   If skipfile is set to 1, sbatch will not attempt to interact with the job
%   file (eg for load balancing). This is useful if the file hasn't been
%   written yet but you still want to queue the job.

% Set some defaults if required
job.host = 'local';
if ~isfield(job,'nodes')
    job.nodes = 1;
end

if ~isfield(job,'ppn')
    job.ppn = 4;
end
if job.ppn>4
    job.ppn = 4;
end

job.ntasks = job.nodes*job.ppn;

if ~isfield(job,'walltime')
    job.walltime = 12.0;
end

% job.mem = 63900*job.nodes;  % Updated 30/08/12 for reduced memory on Darwin
% job.interactive = 0;

% Set up the output filename and log filename
outname = strrep(job.filename,'.hdf5','Output');
outname = strrep(outname,'Input','');
logname = strrep(outname,'Output','');
logname = [logname 'log.txt'];
%logname = ['log' logname '.txt'];

% Preprocess the walltime
hours = floor(job.walltime);
minutes = round((job.walltime-hours)*60);

hours = int2str(hours);
if length(hours) < 2
    hours = ['0' hours];
end

minutes = int2str(minutes);
if length(minutes) < 2
    minutes = ['0' minutes];
end

% Load balance the grid, if necessary

% Choose which service level to run with
if (isfield(job,'SL') && job.SL == 2) || job.walltime > 12.0 || job.nodes > 6
    selector = 'HALL-PULLAN-IFAN-SL2-GPU';%'HALL-SL2-GPU';
elseif  (isfield(job,'SL') && job.SL == 4)
    selector = 'HALL-SL3-GPU';
else
    selector = 'HALL-SL3-GPU';
end

intr = 0;
if job.walltime<=1
    [~,result] = system('squeue -u tsw30 -o "%.5Q"');
    result = strsplit(result,'\n');
    result = result(2:end-1);
    if ~isempty(result)
        pr = str2double(result);
        if any(pr>=9000)
            intr = 0;
        else
            intr = 1;
        end
    else
        intr = 1;
    end
end

% Create the MPI submission script by modifying the template
fidTemplate = fopen('slurmSubmitTemplate_ampere','r');

fidOutput = fopen([job.directory 'slurmSubmitOutput'],'w+');
while ~feof(fidTemplate)
    row = fgetl(fidTemplate);
    row = strrep(row, 'ACCOUNT', selector);
    row = strrep(row, 'NUMBERNODES', int2str(job.nodes));
    row = strrep(row, 'NUMBERTASKS', int2str(job.ntasks));
    row = strrep(row, 'PROCESSES_PER_NODE', int2str(job.ppn));
    row = strrep(row, 'WALLTIME', [hours ':' minutes ':00']);
    row = strrep(row, 'INPUTFILENAME', job.filename);
    row = strrep(row, 'OUTPUTFILENAME', outname);
    row = strrep(row, 'WORKDIRECTORY', job.directory);
    row = strrep(row, 'LOGFILENAME', logname);
    row = strrep(row, 'gpujob', job.name);
    if intr==1
        row = strrep(row, '#!SBATCH --requeue', '#SBATCH --qos=INTR');
    end
    fprintf(fidOutput,'%s\n',row);
end

% Submit to the queue
[~, result] = system(['sbatch ' job.directory 'slurmSubmitOutput']);
% [~,~] = fprintf('exit');
pattern = '\d{3,12}';
fprintf(['\n' result '\n'])

job.id = regexp(result,pattern,'match');
job.id = str2double(job.id{1});

job.logname = logname;
job.inname = job.filename;
job.outname = [outname '.hdf5'];
job.avgname = strrep(job.outname,'.hdf5','_avg.hdf5');
job.remote = 'local';
job.user = 'tsw30';
if ~isfield(job,'localdir')
    job.localdir = job.directory;
end
job.started = 0;
job.finished = 0;

fclose(fidTemplate);


end