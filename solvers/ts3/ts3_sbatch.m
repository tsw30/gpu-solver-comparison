function job = ts3_sbatch(job)
%ts3_sbatch Write a SLURM submission script from a template and submit it.
%   Required job fields: name, filename, directory.
%   Optional fields with fixed-cluster defaults: nodes, ppn, walltime, user, localdir.

%% Defaults
job.host   = 'local';
job.remote = 'local';
if ~isfield(job,'nodes'),    job.nodes    = 1;      end
if ~isfield(job,'ppn'),      job.ppn      = 4;      end
job.ppn     = min(job.ppn,4);                       % max 4 GPUs per node
job.ntasks  = job.nodes * job.ppn;
if ~isfield(job,'walltime'), job.walltime = 12.0;   end
if ~isfield(job,'user'),     job.user     = getenv('USER'); end
if ~isfield(job,'localdir'), job.localdir = job.directory; end

%% Output and log file names
% 'InputFoo.hdf5' -> outname 'FooOutput', logname 'Foolog.txt'
outname = strrep(job.filename,'.hdf5','Output');
outname = strrep(outname,'Input','');
logname = [strrep(outname,'Output','') 'log.txt'];

%% Walltime as HH:MM:SS
hh = floor(job.walltime);
mm = round((job.walltime - hh)*60);
walltime = sprintf('%02d:%02d:00',hh,mm);

%% Interactive queue: use qos=INTR for short single-node jobs,
% provided no other job of ours is queued/running with priority > 9000
intr = job.nodes <= 1 && job.walltime <= 1.0;
if intr
    [~, r] = system(['squeue -u ' job.user ' -h -o "%Q"']);
    if any(str2double(strsplit(strtrim(r))) > 9000)
        intr = false;
    end
end

%% Build submission script from template
fi = fopen('slurmSubmitTemplate_ampere','r');
if fi == -1
    error('ts3_sbatch: could not open template slurmSubmitTemplate_ampere');
end
fo = fopen([job.directory 'slurmSubmitOutput'],'w+');
while ~feof(fi)
    row = fgetl(fi);
    row = strrep(row,'ACCOUNT',            'HALL-AZTEC-SL2-GPU');
    row = strrep(row,'NUMBERNODES',        int2str(job.nodes));
    row = strrep(row,'NUMBERTASKS',        int2str(job.ntasks));
    row = strrep(row,'PROCESSES_PER_NODE', int2str(job.ppn));
    row = strrep(row,'WALLTIME',           walltime);
    row = strrep(row,'INPUTFILENAME',      job.filename);
    row = strrep(row,'OUTPUTFILENAME',     outname);
    row = strrep(row,'WORKDIRECTORY',      job.directory);
    row = strrep(row,'LOGFILENAME',        logname);
    row = strrep(row,'gpujob',             job.name);
    if intr
        row = strrep(row,'#SBATCH --requeue','#SBATCH --qos=INTR');
    end
    fprintf(fo,'%s\n',row);
end
fclose(fo);
fclose(fi);

%% Submit
[~, result] = system(['sbatch ' job.directory 'slurmSubmitOutput']);
fprintf('\n%s\n',result);
job.id = str2double(regexp(result,'\d{3,12}','match','once'));

%% Result file names
job.logname = logname;
job.inname  = job.filename;
job.outname = [outname '.hdf5'];
job.avgname = strrep(job.outname,'.hdf5','_avg.hdf5');

%% Status
job.started  = 0;
job.finished = 0;
end