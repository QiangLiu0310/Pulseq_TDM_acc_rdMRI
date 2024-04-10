function varargout = read_meas_dat(wildfile, user_options, iRep, NReps)
%READ_MEAS_DAT  read in Siemens format raw data, VB-, VD-. and VE-style "meas.dat"
%
% meas_struct = read_meas_dat(filename, <options>)
% meas_struct = read_meas_dat(filename, <options>, <iRep>, <NReps>)
% options_struct = read_meas_dat()
%
%   N.B. "filename" can contain wildcards, such as: 'meas*.dat'
%
%
% "options" is a *structure* which allows changing of some of the defaults:
%
%    options.ApplyFFTScaleFactors         -- set to 0 or 1  (default is 0)
%    options.CanonicalReorderCoilChannels -- set to 0 or 1  (default is 0)
%    options.DispProgress                 -- set to 0 or 1  (default is 1)
%    options.ExtractRepetition            -- set to 0 or 1  (default is 0)
%    options.FlipReflectedLines           -- set to 0 or 1  (default is 1)
%    options.MatrixDouble                 -- set to 0 or 1  (default is 0)
%    options.OfflineSave                  -- set to 0 or 1  (default is 1)
%    options.PhascorCollapseSegments      -- set to 0 or 1  (default is 0)
%    options.ReadMultipleRepetitions      -- set to 0 or 1  (default is 1)
%    options.ReturnCellArray              -- set to 0 - 16  (default is 0)  [REP == 7, LIN == 2]
%    options.ReturnPATRefOnly             -- set to 0 or 1  (default is 0)
%    options.ReturnPhascorOnly            -- set to 0 or 1  (default is 0)
%    options.ReturnStruct                 -- set to 0 or 1  (default is 1)
%    options.SkipRTFeedback               -- set to 0 or 1  (default is 0)
%    options.SMSRefscan                   -- set to 0 or 1  (default is 0)
%    options.SqueezeChannels              -- set to 0 or 1  (default is 1)  [formerly "RemapChannels"]
%    options.StoreTimeStamp               -- set to 0 or 1  (default is 0)
%
% if a field does not exist, then the default is used. "options" is optional :).
%
%
% the raw k-space data is returned in a 16-dimensional array, following
% the convention used in ICE. the mapping between the array dimensions
% and the loopcounters used in ICE is as follows:
%
%    #01:  ColMeas
%    #02:  LinMeas
%    #03:  ChaMeas
%    #04:  SetMeas
%    #05:  EcoMeas
%
%    #06:  PhsMeas
%    #07:  RepMeas
%    #08:  SegMeas
%    #09:  ParMeas
%    #10:  SlcMeas
%
%    #11:  IdaMeas
%    #12:  IdbMeas
%    #13:  IdcMeas
%    #14:  IddMeas
%    #15:  IdeMeas
%
%    #16:  AveMeas
%
%
% ---- TO READ ARBITRARY REPETITIONS: ----
%  - set "options.ExtractRepetition" to 1
%  - specify inputs "iRep" and "NReps" to choose repetitions
%
% note that when "options.ExtractRepetition" is enabled auxiliary data
% stored in the meas.dat file (e.g., PAT refscans) will not be returned. the
% recommended procedure is to (a) read the auxiliary data and first
% repetition by DISABLING "options.ReadMultipleRepetitions", then (b) read
% later repetitions by setting "options.ExtractRepetition" and iRep = 2...N.
%

% (based on the FAMOUS "read_mdh_adc.m" by anders, then mukund, then andre.)

% 2006/dec/04: added support for EPI data
%                - (PHASCOR) extract phase correction lines
%                - (REFLECT) reorder samples from reflected lines

% 2006/dec/06: added support for iPAT data
%                - (NOISEADJSCAN)
%                - (PATREFSCAN)

% 2007/jan/01: added support for 3D EPI phase correction lines

% 2007/mar/29: added canonical preamp-based coil ordering (cf. ICE code)

% 2007/mar/30: implemented Anastasia's suggestions
%                - (PATREFANDIMASCAN) for non-EPI iPAT data

% 2007/may/02: added support for phase stabilization scans
%                - (REFPHASESTABSCAN)
%                - (PHASESTABSCAN)

% 2007/jul/08: added fix for non-contiguous coil channels found in 7T data

% 2007/jul/09: improved support for "meas.out" files lacking headers

% 2007/aug/07: added ability to return all data as fields of a struct

% 2007/aug/14: fixed phase stabilization support for multiecho acquisitions

% 2007/oct/16: fixed order swapping for even and odd phascor lines in EPI

% 2007/oct/19: added ability to return data even if file incomplete

% 2011/oct/24: added support for VD11 data (AY)

% 2011/nov/20: added support for SMS reference data

% 2011/dec/05: added support for returning OFFLINE scans

% 2012/jan/21: added ReturnCellArray option for large, sparse data sets
%
%  to convert the meas.data field from a cell array into a full array,
%  knowledge of exactly how the data is stored in the cell array will be
%  needed. for example, if the acquisition was accelerated in the phase
%  encoding direction such that only the lines specified in "datlines1" were
%  acquired, the conversion through cell2mat() would be
%
%    datprune = cell2mat(meas.data(:,datlines1,:,:,:,:,:,:,:,:,:,:));
%
%  in the special case where no lines are skipped and all cells are
%  populated, the conversion will only require a call to cell2mat(), i.e.,
%
%    dat = cell2mat(meas.data);

% 2012/feb/20: added support for VD11 "multi-RAID" files

% 2012/mar/13: added StoreTimeStamp feature


% Copyright © 2006-2013 Jonathan R. Polimeni and
%   The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense


% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 10/04/2006
% $Id: read_meas_dat.m,v 1.27 2015/07/24 03:19:33 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.27 $';
  if ( nargin == 0 && nargout == 0 ), help(mfilename); return; end;


  %------------------------------------------------------------------------%
  % check version

  if ( mrir_sysutil__versioncheck('7.2') ),
    disp(sprintf('"%s" only supported for MATLAB 7.2 and higher', mfilename));
    return;
  end;


  %------------------------------------------------------------------------%
  % parse options and set defaults

  % the pre-release version installed on the VB13 128-channel host differs
  % slightly from the full-release version, and files must be read in
  % differently depending on the version. yuck!
  IS__VB13_PRE_RELEASE_VERSION = 0;
  IS__RELEASE_VERSION = 1;
  IS__VD_PLATFORM = 0;
  IS__DIMENSION_INFO_PRESENT = 1;

  IS__MULTIRAID = 0;

  DO__HDR_SAVE = 0;

  DO__RECOVER_FROM_INCOMPLETE = 1;

  DO__PREALLOCATE = 0;


  DO__APPLY_FFT_SCALEFACTORS          = 0;
  DO__CANONICAL_REORDER_COIL_CHANNELS = 0;
  DO__EXTRACT_REPETITION              = 0;
  DO__FLIP_REFLECTED_LINES            = 1;
  DO__MDH_SAVE                        = 0;
  DO__PHASCOR_COLLAPSE_SEGMENTS       = 0;
  DO__READ_MULTIPLE_REPETITIONS       = 1;
  DO__RETURN_DEPENDENCY_MEAS          = 0;
  DO__SELECT_CHANNELS                 = 0;
  DO__SQUEEZE_CHANNELS                = 1;

  DO__STORE_TIMESTAMP                 = 0;

  % for MGH version of GESTE sequence -- TODO: automatically detect TENC ACS data
  DO__TENC_ACS                        = 0;

  VERBOSE                             = 0;
  DO__DISP_PROGRESS                   = 1;

  % by default, return individual arrays (for backward compatibility)
  DO__RETURN_STRUCT                   = 1;

  % by default, return all data as full arrays of type SINGLE
  DO__MATRIX_DOUBLE                   = 0;
  DO__MATRIX_SPARSE                   = 0;
  DO__RETURN_CELL_ARRAY               = 0;
  DO__RETURN_PATREF_ONLY              = 0;
  DO__RETURN_PHASCOR_ONLY             = 0;

  % option to return lines that are missing the "ONLINE" flag needed for ICE
  % processing
  DO__SAVE_OFFLINE                    = 1;

  % some sequences, like EPI with DORK, have navigators with
  % RTFeedback flag set and these lines should be skipped, whereas other
  % sequences, like mocoMPRAGE, have all lines flagged as RTFeedback and
  % these lines should be kept.
  DO__SKIP_RTFEEDBACK                 = 0;

  % mechanism for user to specify whether this data contains Simultaneous
  % MultiSlice reference data (since currently there is no flags or
  % conventions to distinguish these lines from other lines!)
  DO__SMSREFSCAN                      = 0;


  % if the "options" struct is provided by caller, then override default
  % flags

  if ( exist('user_options', 'var') ),

    if ( isfield(user_options, 'ApplyFFTScaleFactors') ),
      DO__APPLY_FFT_SCALEFACTORS = user_options.ApplyFFTScaleFactors;
      disp_optional(sprintf(' :APPLY_FFT_SCALEFACTORS = %d', DO__APPLY_FFT_SCALEFACTORS));
    end;

    if ( isfield(user_options, 'CanonicalReorderCoilChannels') ),
      DO__CANONICAL_REORDER_COIL_CHANNELS = user_options.CanonicalReorderCoilChannels;
      disp_optional(sprintf(' :CANONICAL_REORDER_COIL_CHANNELS = %d', DO__CANONICAL_REORDER_COIL_CHANNELS));
    end;

    if ( isfield(user_options, 'DispProgress') ),
      DO__DISP_PROGRESS = user_options.DispProgress;
      disp_optional(sprintf(' :DISP_PROGRESS = %d', DO__DISP_PROGRESS));
    end;

    if ( isfield(user_options, 'ExtractRepetition') ),
      DO__EXTRACT_REPETITION = user_options.ExtractRepetition;
      disp_optional(sprintf(' :EXTRACT_REPETITION = %d', DO__EXTRACT_REPETITION));
      if ( nargin < 3 ),
        iRep = 2;
      end;
      if ( nargin < 4 ),
        NReps = 1;
      end;
      if ( iRep < 2 ),
        disp(sprintf('==> [%s]: to read 1st rep., set ReadMultipleRepetitions = 0 -- see help ("TO READ ARBITRARY REPETITIONS")', mfilename))
      end;
    end;

    if ( isfield(user_options, 'FlipReflectedLines') ),
      DO__FLIP_REFLECTED_LINES = user_options.FlipReflectedLines;
      disp_optional(sprintf(' :FLIP_REFLECTED_LINES = %d', DO__FLIP_REFLECTED_LINES));
    end;

    if ( isfield(user_options, 'MatrixDouble') ),
      DO__MATRIX_DOUBLE = user_options.MatrixDouble;
      disp_optional(sprintf(' :MATRIX_DOUBLE = %d', DO__MATRIX_DOUBLE));
    end;

    if ( isfield(user_options, 'OfflineSave') ),
      DO__SAVE_OFFLINE = user_options.OfflineSave;
      disp_optional(sprintf(' :SAVE_OFFLINE = %d', DO__SAVE_OFFLINE));
    end;

    if ( isfield(user_options, 'MDHSave') ),
      DO__MDH_SAVE = user_options.MDHSave;
      disp_optional(sprintf(' :MDH_SAVE = %d', DO__MDH_SAVE));
    end;

    if ( isfield(user_options, 'PhascorCollapseSegments') ),
      DO__PHASCOR_COLLAPSE_SEGMENTS = user_options.PhascorCollapseSegments;
      disp_optional(sprintf(' :PHASCOR_COLLAPSE_SEGMENTS = %d', DO__PHASCOR_COLLAPSE_SEGMENTS));
    end;

    if ( isfield(user_options, 'ReadMultipleRepetitions') ),
      DO__READ_MULTIPLE_REPETITIONS = user_options.ReadMultipleRepetitions;
      disp_optional(sprintf(' :READ_MULTIPLE_REPETITIONS = %d', DO__READ_MULTIPLE_REPETITIONS));
    end;

    if ( isfield(user_options, 'ReturnCellArray') ),
      DO__RETURN_CELL_ARRAY = user_options.ReturnCellArray;
      disp_optional(sprintf(' :RETURN_CELL_ARRAY = %d', DO__RETURN_CELL_ARRAY));
    end;

    if ( isfield(user_options, 'ReturnPATRefOnly') ),
      DO__RETURN_PATREF_ONLY = user_options.ReturnPATRefOnly;
      disp_optional(sprintf(' :RETURN_PATREF_ONLY = %d', DO__RETURN_PATREF_ONLY));
    end;

    if ( isfield(user_options, 'ReturnPhascorOnly') ),
      DO__RETURN_PHASCOR_ONLY = user_options.ReturnPhascorOnly;
      disp_optional(sprintf(' :RETURN_PHASCOR_ONLY = %d', DO__RETURN_PHASCOR_ONLY));
    end;

    if ( isfield(user_options, 'ReturnDependencyMeas') ),
      DO__RETURN_DEPENDENCY_MEAS = user_options.ReturnDependencyMeas;
      disp_optional(sprintf(' :RETURN_DEPENDENCY_MEAS = %d', DO__RETURN_DEPENDENCY_MEAS));
    end;

    if ( isfield(user_options, 'ReturnStruct') ),
      DO__RETURN_STRUCT = user_options.ReturnStruct;
      disp_optional(sprintf(' :RETURN_STRUCT = %d', DO__RETURN_STRUCT));
    end;

    if ( isfield(user_options, 'SelectChannels') ),
      if ( ~isempty(user_options.SelectChannels) ),
	DO__SELECT_CHANNELS = 1;
	coil_selected = user_options.SelectChannels;
	disp_optional(sprintf(' :SELECT_CHANNELS = %d', DO__SELECT_CHANNELS));
      end;
    end;
    
    if ( isfield(user_options, 'SkipRTFeedback') ),
      DO__SKIP_RTFEEDBACK = user_options.SkipRTFeedback;
      disp_optional(sprintf(' :SKIP_RTFEEDBACK = %d', DO__SKIP_RTFEEDBACK));
    end;

    if ( isfield(user_options, 'SMSRefscan') ),
      DO__SMSREFSCAN = user_options.SMSRefscan;
      disp_optional(sprintf(' :SMSREFSCAN = %d', DO__SMSREFSCAN));
    end;

    if ( isfield(user_options, 'SqueezeChannels') ),
      DO__SQUEEZE_CHANNELS = user_options.SqueezeChannels;
      disp_optional(sprintf(' :SQUEEZE_CHANNELS = %d', DO__SQUEEZE_CHANNELS));
    end;

    if ( isfield(user_options, 'RemapChannels') ),
      disp('');
      disp('------------------------------------------------------------------------------');
      warning(' !!! option "RemapChannels" has been renamed to "SqueezeChannels" !!! ');
      disp('------------------------------------------------------------------------------');
      disp('');
    end;

    if ( isfield(user_options, 'StoreTimeStamp') ),
      DO__STORE_TIMESTAMP = user_options.StoreTimeStamp;
      disp_optional(sprintf(' :STORE_TIMESTAMP = %d', DO__STORE_TIMESTAMP));
    end;

    if ( isfield(user_options, 'TemporallyEncodedACS') ),
      DO__TENC_ACS = user_options.TemporallyEncodedACS;
      disp_optional(sprintf(' :TENC_ACS = %d', DO__TENC_ACS));
    end;

    if ( isfield(user_options, 'VERBOSE') ),
      VERBOSE = user_options.VERBOSE;
      disp_optional(sprintf(' :VERBOSE = %d', VERBOSE));
    end;

    user_option_list = fieldnames(user_options);

  else,
    user_option_list = {};
  end;

  options = struct('ApplyFFTScaleFactors',          DO__APPLY_FFT_SCALEFACTORS, ...
                   'CanonicalReorderCoilChannels',  DO__CANONICAL_REORDER_COIL_CHANNELS, ...
                   'DispProgress',                  DO__DISP_PROGRESS, ...
                   'ExtractRepetition',             DO__EXTRACT_REPETITION, ...
                   'FlipReflectedLines',            DO__FLIP_REFLECTED_LINES, ...
                   'MatrixDouble',                  DO__MATRIX_DOUBLE, ...
                   'MDHSave',                       DO__MDH_SAVE, ...
                   'OfflineSave',                   DO__SAVE_OFFLINE, ...
                   'PhascorCollapseSegments',       DO__PHASCOR_COLLAPSE_SEGMENTS, ...
                   'ReadMultipleRepetitions',       DO__READ_MULTIPLE_REPETITIONS, ...
                   'ReturnCellArray',               DO__RETURN_CELL_ARRAY, ...
                   'ReturnPATRefOnly',              DO__RETURN_PATREF_ONLY, ...
                   'ReturnPhascorOnly',             DO__RETURN_PHASCOR_ONLY, ...
                   'ReturnDependencyMeas',          DO__RETURN_DEPENDENCY_MEAS, ...
                   'ReturnStruct',                  DO__RETURN_STRUCT, ...
                   'SelectChannels',                DO__SELECT_CHANNELS, ...
                   'SkipRTFeedback',                DO__SKIP_RTFEEDBACK, ...
                   'SMSRefscan',                    DO__SMSREFSCAN, ...
                   'SqueezeChannels',               DO__SQUEEZE_CHANNELS, ...
                   'StoreTimeStamp',                DO__STORE_TIMESTAMP, ...
                   'TemporallyEncodedACS',          DO__TENC_ACS, ...
                   'VERBOSE',                       VERBOSE);

  if ( ~isempty(user_option_list) ),

    user_option_list_VALID = isfield(options, user_option_list);

    if ( ~all(user_option_list_VALID) ),
      invalid_items = find(~user_option_list_VALID);
      for item = 1:length(invalid_items),
        warning(sprintf('option "%s" not recognized', strvcat(user_option_list{invalid_items(item)})));
      end;
    end;

  end;

  % feature: if user does not specify a file but asks for a return argument,
  % return the default options struct. this might be helpful for the user to
  % see all the options and to avoid having to type the field names
  % manually.
  if ( nargin == 0 && nargout == 1 ),
    disp(sprintf('==> [%s]: returning default options as struct', mfilename));
    varargout{1} = options;
    return;
  end;


  %------------------------------------------------------------------------%
  % ICE constants

  % from "MrServers/MrVista/include/Ice/IceDefs.h":
  ICE_RAWDATA_SCALE       = 131072.0;  % 64 ^ 3 / 2
  K_ICE_AMPL_SCALE_FACTOR = 80 * 20 * ICE_RAWDATA_SCALE / 65536;


  %------------------------------------------------------------------------%
  % check for file

  % if no input, then simply return options struct
  if ( isempty(wildfile) ),
    disp('no input file detected -- returning "options" struct only');
    varargout{1} = options;
    return;
  end;

  global DATADIR
  if ( ~isempty(DATADIR) ),
    if ( ~exist(DATADIR, 'dir') ), DATADIR = ''; end;
  end;

  if ( iscellstr(wildfile) ),
    wildfile = cell2mat(wildfile);
  end;

  % expand wild cards
  pathstr = fileparts(wildfile);

  if ( isempty(pathstr) && ~isempty(DATADIR) ),
    disp(sprintf('(setting path to DATADIR "%s")', DATADIR));
    pathstr = DATADIR;
    wildfile = fullfile(pathstr, wildfile);
  end;

  if ( isempty(dir(wildfile)) ),
    error('pattern [%s] does not match existing meas file\n', wildfile);
  end;

  matching_files = dir(wildfile);

  % in case multiple files match, extract first file
  measfile = fullfile(pathstr, getfield(matching_files, {1}, 'name'));

  % report matching file name if regexp found
  if ( ~isempty(regexp(wildfile, '*')) || ~isempty(regexp(wildfile, '?')) ),
    disp(sprintf('reading file "%s"...', measfile));
  end;

  if ( ~exist(measfile, 'file') ),
    error('file [%s] does not exist', measfile);
  end;


  %========================================================================%
  % B E G I N    R E A D I N G
  %========================================================================%

  t0 = clock;

  [pathstr, filestr, extstr] = fileparts(measfile);

  [fp, errstr] = fopen(measfile, 'r', 'l');
  if ( fp == -1 ),
    error(errstr);
  end;


  % determine size (in bytes) of ascii header files stored in the 'meas.dat'
  % format (i.e., "Config_.evp", "Dicom_.evp", etc) to skip over them all.
  % [note: this works for VB11A also---the first integer is 32, which is
  % the number of bytes to skip at the beginning!]
  data_start = fread(fp, 1, 'uint32');

  if ( data_start == 0 ),    % Is this a VD/VE data file?
    IS__VD_PLATFORM = 1;
  end;

  if ( IS__VD_PLATFORM ),

    % count number of concatenated meas.dat files stored in this file
    num_dependent_files = fread(fp, 1, 'uint32');

    if ( num_dependent_files > 1 ),
      IS__MULTIRAID = 1;
      fprintf(1, '\n ***MULTI-RAID file detected***\n  %d meas contained within file \n\n', num_dependent_files);

      if ( DO__RETURN_DEPENDENCY_MEAS == 0 ),
        fprintf(1, ' input "%s" contains [%d] dependent meas.dat ADJUSTMENT file(s)\n', strcat(filestr, extstr), num_dependent_files-1);
        disp(' reading in *main data file* in multi-RAID file...');
        if ( num_dependent_files == 2 ),
          fprintf(1, ' re-run and set option field "ReturnDependencyMeas" to 1 to load another meas within this file\n\n');
        elseif ( num_dependent_files == 3 ),
          fprintf(1, ' re-run and set option field "ReturnDependencyMeas" to 1 or 2 to load another meas within this file\n\n');
        else,
          fprintf(1, ' re-run and set option field "ReturnDependencyMeas" to 1--%d to load another meas within this file\n\n', num_dependent_files-1);
        end;
      else,
        fprintf(1, ' reading meas %d of %d within the multi-RAID file\n\n', options.ReturnDependencyMeas, num_dependent_files);
      end;
    end;

    if ( ~IS__MULTIRAID ),
      meas_file_within_multifile = 0;
    end;

    % by default read the last file, which is the main file stored after the adjustment files
    if ( DO__RETURN_DEPENDENCY_MEAS == 0 ),
      meas_file_within_multifile = (num_dependent_files - 1);
    else,
      meas_file_within_multifile = options.ReturnDependencyMeas - 1;
    end;

    % magic number alert: 152
    fseek(fp, 16+(meas_file_within_multifile*152), 'bof');
    header_start = fread(fp, 1, 'uint64');
    fseek(fp, header_start, 'bof');
    data_start = header_start + fread(fp, 1, 'uint32');
  else,
    header_start = 0;
  end;


  % if file has no header (e.g., if its in "meas.out" format),
  % "data_start" should be 32.
  if ( data_start == 32 ),
    disp_optional('no header detected! assuming "meas.out" format...');

    IS__RELEASE_VERSION = 0;

    [pathstr, name, ext] = fileparts(measfile);
    meas_asc = fullfile(pathstr, [name, '.asc']);

    if ( exist(meas_asc, 'file') ),
      disp_optional('discovered corresponding "meas.asc" file! parsing...');

      [asc, errstr] = fopen(meas_asc, 'r', 'l');
      if ( asc == -1 ), error(errstr); end;

      % read header into one string for parsing
      header = fscanf(asc, '%c');

      % HACK: sometimes two headers are stored in one file (e.g.,
      % AdjCoilSens); find the end the YAPS marker and if there are two
      % truncate header at the first one
      header_end = regexp(header, '}\n### ASCCONV BEGIN');  % avoids MeasYaps
      if ( length(header_end) > 1 ),
        header = header(1:header_end(end));
      end;

      fclose(asc);

    else,
      % set header to empty string
      header = '';
    end;

    % can't sort channels without header information  :(
    DO__CANONICAL_REORDER_COIL_CHANNELS = 0;

    % jump to beginning of binary data, let's get to work!
    fseek(fp, data_start, 'bof');

  elseif (  ( DO__PREALLOCATE || DO__CANONICAL_REORDER_COIL_CHANNELS || DO__RETURN_STRUCT )  ),

    % read header into one string for parsing
    header = char(fread(fp, data_start-header_start-4, 'uchar').');

    param_list = {'NColMeas', 'NLinMeas', 'NChaMeas', 'NSetMeas', 'NEcoMeas', ...
                  'NPhsMeas', 'NRepMeas', 'NSegMeas', 'NParMeas', 'NSlcMeas', ...
                  'NIdaMeas', 'NIdbMeas', 'NIdcMeas', 'NIddMeas', 'NIdeMeas', ...
                  'NAveMeas'};

    dimensions = cell2struct(cell(length(param_list),1), param_list, 1);
    dim = [];

    % scan through header for each of the ICE dimension values
    for ind = 1:length(param_list),
      param = param_list{ind};

      %%% SPECIAL CASE: "NSegMeas"

%%      % the number of segments is listed in two places in the header with the
%%      % field names "NSegMeas" and "NSeg", and for some reason only the "NSeg"
%%      % field gives the correct number of segments; SO for this field we break
%%      % with the convention
%%      % UPDATE: it appears that "NSegMeas" corresponds to the cumulative
%%      % number of distinct segments appearing in the loop counters, whereas
%%      % "NSeg" is the true number of segments that is the same as the
%%      % number of shots, i.e., they should differ by a factor of 2 when the
%%      % "OnlineTSE" functor is in use.
%%      if ( strcmp(param, 'NSegMeas') ),
%%        param = 'NSeg';
%%      end;

      % exploit MATLAB regexp machinery to pull out parameter/value pairs
      match = regexp(header, ['(?<param>' param, ').{0,5}\{\s*(?<value>\d*)\s*\}'], 'names');

      % check if no match is found
      if ( isempty(match) ),
        if ( IS__DIMENSION_INFO_PRESENT ),
          IS__DIMENSION_INFO_PRESENT = 0;
          warning('SIEMENS:IO:versioning', 'missing data dimension info in header');
        end;
        continue;
      end;

      % consider only last match (there can be as many as three in Config_.evp)
      match = match(end);

      % empty means number of elements in this dimension = 1
      if ( isempty(match.value) ),
        match.value = '1';
      end;

      % save out struct and numerical array
      dim(ind) = str2double(match.value);
      dimensions.(param_list{ind}) = dim(ind);

    end;

    % VB11A workaround hack (for EPI only)
    if ( ~IS__RELEASE_VERSION ),
      dimensions.NColMeas = 0;
      dimensions.NAveMeas = 1;
      dimensions.NSegMeas = 2;
    end;


    %------------------------------------------------------------------------%
    % extract FFT scalefactors (if utility file is in path)

    if ( DO__APPLY_FFT_SCALEFACTORS && DO__EXTRACT_REPETITION ),
      DO__APPLY_FFT_SCALEFACTORS = 0;
      disp_optional('cannot apply FFT scale factors when extracting a repetition due to lack of headers -- skipping scaling');
    end;

    if ( DO__APPLY_FFT_SCALEFACTORS && ...
         exist('read_meas_dat__fft_scalefactors', 'file') ),
      [fft_scale, fft_scale_channel] = read_meas_dat__fft_scalefactors(header);
      if ( ~isempty( fft_scale ) ),
        disp_optional('applying FFT scale factors.');
      else,
        DO__APPLY_FFT_SCALEFACTORS = 0;
        disp_optional('FFT scale factors not found!!!');
      end;
    else,
      fft_scale = [];
      disp_optional('ignoring FFT scale factors.');
    end;


    %------------------------------------------------------------------------%
    % compute canonical coil ordering from coil element strings---the channel
    % number assigned to each coil is determined by the order in which the
    % channels are selected in SYNGO, and NOT by the preamp channel numbers,
    % so to enforce a consistent ordering across scans that is independent of
    % the operator's coil selections we can sort by the fixed preamp channel
    % strings.

    % NOTE: this issue tends to arise only for customer coil files; from the
    % cases we've seen product coil files tend to be more sane and (a)
    % select groups of elements with a single button and (b) order the
    % channels consistently. but cannot prove the null hypothesis---maybe
    % product coil files have this issue sometimes too.

    if ( DO__CANONICAL_REORDER_COIL_CHANNELS && DO__EXTRACT_REPETITION ),
      DO__CANONICAL_REORDER_COIL_CHANNELS = 0;
      disp_optional('cannot reorder channels when extracting a repetition due to lack of headers -- skipping reordering');
    end;

    if ( DO__CANONICAL_REORDER_COIL_CHANNELS && ...
         exist('read_meas_dat__reorder_coil_channels', 'file') ),
      disp_optional('reordering coil channels to canonical preamp-based ordering.');

      [coil_index, coil_order] = read_meas_dat__reorder_coil_channels(header);

      % to map using index:
      %   dataval(coil_index);
      % to map using order:
      %   coil_order(dataind);

      if ( DO__APPLY_FFT_SCALEFACTORS ),
        if ( length(fft_scale) ~= length(coil_index) ),
          DO__APPLY_FFT_SCALEFACTORS = 0;
          disp_optional('mismatch between number of channels and FFT scale factors!!!  -->  ignoring FFT scale factors...');
        else,
          % don't forget to reorder FFT scalefactors!
          fft_scale = fft_scale(coil_index);
        end;
      end;

    else,
      % clear flag if auxiliary function not in path
      DO__CANONICAL_REORDER_COIL_CHANNELS = 0;
    end;

    %------------------------------------------------------------------------%

  end;  % end of header parsing when present, if..else..end

  % jump past header to beginning of binary data
  fseek(fp, data_start, 'bof');

  % save current file position, then jump to end to calculate remaining data size
  fpos = ftell(fp);
  fseek(fp, 0, 'eof');
  eof = ftell(fp);
  databytes = eof;

  % return to saved file position
  fseek(fp, fpos, 'bof');

  disp_optional(sprintf('binary data after header: %10.2f MB', databytes/2^20));


  % initializations...
  meas_num = 0;
  scan_num = 0;

  maxLin = -1;
  maxRep = -1;

  %mdh = read_meas__mdh_struct;

  %--------------------------------------------------------------------------%

  matrixtype = @single;
  matrixtypestr = 'single';

  % unfortunately, matlab does not support sparse arrays whose number of
  % dimensions is larger than 2.   :(
  if ( DO__MATRIX_SPARSE ),
    warning('sparse matrices are currently unsupported');
    %%%matrixtype = @sparse;
  end;

  if ( DO__MATRIX_DOUBLE ),
    matrixtype = @double;
    matrixtypestr = 'double';
  end;



  meas = struct;
  meas.file = strcat(filestr, extstr);
  meas.readtime = '';

  meas.acqtime = [];

  % number of measurements along all 16 dimensions can be found in the
  % 'Config_.evp' file contained within "meas.dat" header; allocate memory
  % for measurement data, and grow data matrix with each line (this method
  % is, surprisingly, faster than pre-allocation in some cases depending on
  % array shape!)

  if ( DO__PREALLOCATE ),
    tic;
    fprintf(' allocating memory...');
    meas.data              = complex(zeros(struct2array(dimensions), matrixtypestr));
    fprintf('done!  ');
    toc;
    fprintf('\n');
  end;

  if ( DO__RETURN_CELL_ARRAY ),
    meas.data                = cell(1);
    %  meas.data                = {[]};
    patrefscan               = cell(1);
  else,
    meas.data                = matrixtype([]);
    patrefscan               = matrixtype([]);
    %  patrefscan               = {matrixtype([])};
    %  patrefscan               = {[]};
    %  patrefscan               = cell(1);
  end;

  FLAG__data               = 0;
  FLAG__data_reflect       = 0;
  FLAG__data_swapped       = 0;

  FLAG__allow_asym         = 0;

  FLAG__phascor1d          = 0;  % indicating presence of any type of phascor line

  data_phascor1d           = matrixtype([]);
  FLAG__data_phascor1d     = 0;
  data_phascor2d           = matrixtype([]);
%  scan_phascor2d           = matrixtype([]);
  FLAG__data_phascor2d     = 0;
  data_fieldmap            = matrixtype([]);
  FLAG__data_fieldmap      = 0;
  offline                  = matrixtype([]);
  FLAG__offline            = 0;
  noiseadjscan             = matrixtype([]);
  FLAG__noiseadjscan       = 0;
%  scan_patrefscan          = matrixtype([]);
  FLAG__patrefscan         = 0;
  patrefscan_phascor       = matrixtype([]);
  FLAG__patrefscan_phascor = 0;
  FLAG__patrefandimascan   = 0;
  pat_acq_firstscan        = [];

  smsrefscan                  = matrixtype([]);
  FLAG__smsrefscan            = 0;
  smsrefscan_phascor          = matrixtype([]);
  FLAG__smsrefscan_phascor    = 0;

  phasestabtime            = double([]);
  phasestabscan            = matrixtype([]);
  refphasestabtime         = double([]);
  refphasestabscan         = matrixtype([]);
  patrefphasestabscan      = matrixtype([]);
  patrefphasestabtime      = double([]);
  FLAG__refphasestabscan   = 0;
  FLAG__phasestabscan      = 0;


  FLAG__data_phascor1d_orderswap = 0;

  FLAG__patrefscan_phascor_orderswap = 0;
  FLAG__smsrefscan_phascor_orderswap = 0;

  FLAG__repeated_patrefscan = 0;
  FLAG__repeated_smsrefscan = 0;

  FLAG__smsrefscan_completed = 0;


  FLAG__syncdata = 0;

  FLAG__rtfeedback = 0;
  FLAG__dorklines = 0;

 % set to zero if trouble
  FLAG__status_OK = 1;


  timestamp_data = [];
  timestamp_data_phascor1d = [];


  % convert binary mask into strings
  EvalInfoMask = read_meas__evaluation_info_mask_definition;
  EvalInfoNames = fieldnames(EvalInfoMask);

  meas_timestart = 0;

  % channel index offset (needed only if IS__VB13_PRE_RELEASE_VERSION)
  channel_offset = NaN;

  % store all channel indices in case non-contiguous (e.g., 7T host)
  FLAG__stored_channel_indices = 0;
  channel_indices = [];

  FLAG__channel_squeeze = 0;
  channel_index_map = [];

  ACQEND = 0;

  try,
    %------------------------------------------------------------------------%

    ulDMALength = uint16(fread(fp, 1, 'uint16'));
    ulFlags1    = uint8(fread(fp, 1, 'uint8'));
    ulFlags2    = uint8(fread(fp, 1, 'uint8'));
    % overwritten after first MDH is read -- this may not be needed here
    lMeasUID    = int32(fread(fp, 1, 'int32'));

    fseek(fp, fpos, 'bof');

    if ( IS__VD_PLATFORM ),
      mdh_length_float32 = 48;
      mdh_ch_length_float32 = 8;
    else,
      mdh_length_float32 = 32;
    end;

    data_begin = -1;
    firstdata_pos = -1;
    bytes_per_rep = -1;
    FLAG__FOUND_EXTRACT_REPETITION = 0;

    DO__break_from_innerloop = 0;

    %======================================================================%
    % M A I N    L O O P        S T A R T
    %======================================================================%
    % loop through file, peeling off one MDH and ADC line at a time
    while ( ~ACQEND ),

      if ( DO__break_from_innerloop ),
        break;
      end;

      meas_pos = ftell(fp);
      meas_num = meas_num + 1;

      % TODO: save the first MDH for each data type as an example (e.g.,
      % noisescan, patrefscan, data...)
      if ( DO__MDH_SAVE ),
        idx = meas_num;
      else,
        % overwrite mdh(1) with every measurement line
        idx = 1;
      end;

      mdh_binary = fread(fp, mdh_length_float32*4, 'uchar=>uchar');
      [mdh(idx), scan_num] = read_meas_dat__mdh_binary__alt(mdh_binary, 1);

%      if ( mdh(idx).ulTimeSinceLastRF ~= 0 ),
%        mdh(idx).ulTimeSinceLastRF
%        keyboard
%      end;

      aulEvalInfoMask = uint64(mdh(idx).aulEvalInfoMask(1));

      ACQEND = read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_ACQEND);
      if ( ACQEND ),


        % the last MDH and ADC line contains a 16-sample measurement of
        % unknown utility, so its stored in 'adc_cplx' without being
        % overwritten, but is not returned.
        break;

      end;  % IF "ACQEND"


      % for some auxiliary scans (e.g., "AdjCoilSens" for the prescan
      % normalize), the first scan contains gibberish and its mask contains
      % the 'SYNCDATA' bit, so it is likely passing some message back to the
      % sequence that we don't need. if this is found, just skip ahead to
      % the next MDH.
      if ( ~FLAG__syncdata && read_meas__extract_flag(aulEvalInfoMask, 'SYNCDATA') ),
        fread(fp, double(mdh(idx).ulDMALength) - (mdh_length_float32*4), 'uint8');
        meas_num = meas_num - 1;
        continue;
      end;
      FLAG__syncdata = 1;


      % save a copy of the first MDH
      if ( meas_num == 1 ),
        mdh_first = mdh(1);
        meas_timestart = mdh(1).flTimeStamp_ms;
        lMeasUID = mdh(1).lMeasUID;
      end;


      % heuristic test if MDH read was successful
      if ( mdh(idx).lMeasUID ~= lMeasUID ),

        dstr = sprintf('<!> [%s]:  ABORTING read, error detected in meas file', mfilename);
        warning('SIEMENS:IO:format', 'MeasUID mismatch');
        %keyboard

        if ( DO__RECOVER_FROM_INCOMPLETE ),

          fprintf(1, '\n');
          disp_optional(dstr);
          fprintf(1, '\n');
          FLAG__status_OK = 0;

          break;

        else,
          error(dstr);
        end;

      end;


      % make a fix to channel ordering for OLD 128-channel host bug
      if ( IS__VB13_PRE_RELEASE_VERSION ),

        % prior to the full release, a "known bug" in the channel ID
        % numbering conspired to index the channels sequentially BUT began
        % the indexing with a seemingly *random* integer. (twitzel@nmr
        % implemented a similar workaround in his code.)

        FIRSTSCANINSLICE = read_meas__extract_flag(aulEvalInfoMask, 'FIRSTSCANINSLICE');
        if ( FIRSTSCANINSLICE && (mdh(idx).ulScanCounter == 1) ),

          if ( isnan(channel_offset) ),
            channel_offset = mdh(idx).ushChannelId;
          end;

          mdh(idx).ushChannelId = mdh(idx).ushChannelId - channel_offset;

        end;

      end;  % IF "IS__VB13_PRE_RELEASE_VERSION"


      if ( DO__EXTRACT_REPETITION ),

        % skip all lines that are not data lines
        if ( read_meas__extract_bit(aulEvalInfoMask, ...
                                    EvalInfoMask.MDH_NOISEADJSCAN)  || ...
             read_meas__extract_bit(aulEvalInfoMask, ...
                                    EvalInfoMask.MDH_PATREFSCAN)    || ...
             read_meas__extract_bit(aulEvalInfoMask, ...
                                    EvalInfoMask.MDH_PHASESTABSCAN) ),
          if ( IS__VD_PLATFORM ),
            byte_skip__remaining_meas = 0;
            byte_skip__remaining_line = 4*(mdh_ch_length_float32 + double(mdh(idx).ushSamplesInScan)*2) * double(mdh(idx).ushUsedChannels);
          else,
            byte_skip__remaining_meas = 4*[ double(mdh(idx).ushSamplesInScan)*2 ];
            byte_skip__remaining_line = (mdh_length_float32*4 + byte_skip__remaining_meas) * ( double(mdh(idx).ushUsedChannels) - 1 );
          end;

          fseek(fp, byte_skip__remaining_meas + byte_skip__remaining_line, 0);

          continue;
        end;

        % positive once we've seen a line of image data and recorded its position
        if ( firstdata_pos > 0 && ~FLAG__FOUND_EXTRACT_REPETITION ),

          % skip the rest of the first repetition
          if ( mdh(idx).ushRepetition < 1 ),
            if ( IS__VD_PLATFORM ),
              byte_skip__remaining_meas = 0;
              byte_skip__remaining_line = 4*(mdh_ch_length_float32 + double(mdh(idx).ushSamplesInScan)*2) * double(mdh(idx).ushUsedChannels);
            else,
              byte_skip__remaining_meas = 4*[ double(mdh(idx).ushSamplesInScan)*2 ];
              byte_skip__remaining_line = (mdh_length_float32*4 + byte_skip__remaining_meas) * ( double(mdh(idx).ushUsedChannels) - 1 );
            end;

            fseek(fp, byte_skip__remaining_meas + byte_skip__remaining_line, 0);

            continue;
          end;


          % if make it to this spot, we're past first repetition
          firstdata_end_pos = meas_pos;


          % skip all until we find the target repetition
          if ( bytes_per_rep < 0 ),

            fseek(fp, meas_pos, 'bof');

            bytes_per_rep = firstdata_end_pos - firstdata_pos;

            disp_optional(sprintf(' [[ seeking repetition %d... ]]', iRep));
            fseek(fp, bytes_per_rep * (iRep-2), 'cof');

            fpos_rep = ftell(fp);
%      mdh_binary = fread(fp, mdh_length_float32*4, 'uchar=>uchar');
%      [mdh(idx), scan_num] = read_meas_dat__mdh_binary(mdh_binary)
%      fseek(fp, fpos_rep, 'bof');
            channel_indices = [];

            FLAG__FOUND_EXTRACT_REPETITION = 1;

            continue;
          end;

        end;

      end;  % DO__EXTRACT_REPETITION

      if ( DO__EXTRACT_REPETITION && ((mdh(idx).ushRepetition+1) > iRep + NReps - 1) ),
        % all requested repetitions have been found!
        break;
      end;

      % How many channels' data follows this MDH?
      if ( IS__VD_PLATFORM ),           % For VD11 or newer, there is one MDH per line
        nchannel = mdh(idx).ushUsedChannels;
      else,              % For VB17 or older, there is one MDH per line per channel
        nchannel = 1;
      end;

      %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
      for ichannel = 1:nchannel,

        if ( DO__EXTRACT_REPETITION && (firstdata_pos == meas_pos) ),
          continue;
        end;

        if ( IS__VD_PLATFORM ),  % Read channel header
          mdh_binary = fread(fp, mdh_ch_length_float32*4, 'uchar=>uchar');
          mdhch = read_meas_dat__mdh_ch_binary(mdh_binary);

          % heuristic test if MDH read was successful
          if ( mdh(idx).lMeasUID ~= lMeasUID ),

            dstr = sprintf('<!> [%s]:  ABORTING read, error detected in meas file', mfilename);
            warning('SIEMENS:IO:format', 'MeasUID mismatch');

            if ( DO__RECOVER_FROM_INCOMPLETE ),
              fprintf(1, '\n');
              disp_optional(dstr);
              fprintf(1, '\n');
              FLAG__status_OK = 0;

              DO__break_from_innerloop = 1;
              break;

            else,
              error(dstr);
            end;
          end;
          mdh(idx).ushChannelId = mdhch.ushChannelId;
        end;

        % finally, after MDH, read in the data
        adc = (fread(fp, double(mdh(idx).ushSamplesInScan)*2, 'float32=>single'));
        adc_real = adc(1:2:end);
        adc_imag = adc(2:2:end);
        adc_cplx = K_ICE_AMPL_SCALE_FACTOR * complex(adc_real,adc_imag);

        if ( DO__APPLY_FFT_SCALEFACTORS ),
          % scale whatever data comes in (QUESTION: should the noise data be scaled?)
          adc_cplx = adc_cplx * fft_scale( double(mdh(idx).ushChannelId+1) );
        end;

        if ( read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_REFLECT) ),
          if ( ~FLAG__data_reflect ),
            if ( DO__FLIP_REFLECTED_LINES ),
              disp_optional(sprintf(' REFLECT detected (seg #%d), reversing lines.', mdh(idx).ushSeg+1));
            else,
              disp_optional(sprintf(' REFLECT detected (seg #%d).', mdh(idx).ushSeg+1));
            end;
            FLAG__data_reflect = 1;
          end;
          if ( DO__FLIP_REFLECTED_LINES ),
            adc_cplx = flipud(adc_cplx);
          end;
        end;

        if ( read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_RTFEEDBACK) ),
          if ( ~FLAG__rtfeedback ),
            % first occurrence
            FLAG__rtfeedback = 1;
            if ( ~DO__SKIP_RTFEEDBACK ),
              disp_optional(' RTFEEDBACK detected, retaining all feedback lines.');
            else,
              disp_optional(' RTFEEDBACK detected, skipping all feedback lines.');
            end;
          end;  %% first occurrence

          % prior to VD11, all lines with RTFEEDBACK=1 used to be discarded,
          % regardless of the PHASCOR bit. the PHASCOR check below was added
          % for VD11 data, since some Skyra sequences use DORK, which
          % involves phase correction lines to be processed via RTFEEDBACK.
          % these lines should not be discarded.
          if ( read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_PHASCOR) ),
            if ( ~FLAG__dorklines ),
              % first occurrence
              FLAG__dorklines = 1;
              if ( ~DO__SKIP_RTFEEDBACK ),
                disp_optional(' DORK detected.');
              else,
                disp_optional(' DORK detected, now retaining only *PHASCOR* feedback lines.');
              end;
            end; %% first occurrence
          else,
            % if not a DORK line *and* skip option is set, skip this line!
            if ( DO__SKIP_RTFEEDBACK ), continue; end;
          end;
        end;


        if ( read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_SWAPPED) ),
          if ( ~FLAG__data_swapped ),
            % not sure what to do about swapped lines, but for now just report it to user.
            disp_optional(' SWAPPED detected, ignoring.');
            FLAG__data_swapped = 1;
          end;
        end;

        if ( mdh(idx).ushKSpaceCentreColumn ),
          if ( ~FLAG__allow_asym ),

            samples_fullscan = 2*[double(mdh(idx).ushSamplesInScan) - double(mdh(idx).ushKSpaceCentreColumn)];
            samples_skipped = samples_fullscan - double(mdh(idx).ushSamplesInScan);
            samples_asym = double(mdh(idx).ushSamplesInScan) / samples_fullscan;

            if ( samples_asym ~= 1 ),
              disp_optional(sprintf(' ASYMMETRIC ECHO ("allow asymmetry") detected: asym = %2.4f', samples_asym));
              FLAG__allow_asym = 1;
	      ushSamplesInScan = double(mdh(idx).ushSamplesInScan);
            end;

          end;
        end;


        % for testing, abort after first repetition
        if ( ~DO__EXTRACT_REPETITION && ~DO__READ_MULTIPLE_REPETITIONS && mdh(idx).ushRepetition > 0 ...
             && FLAG__data ),

          % [2013/mar/23] the FLAG__data is necessary to accommodate data with multiple ACS repetitions

          %disp_optional('aborting after first repetition...');

          DO__break_from_innerloop = 1;
          break;              % Break from FOR loop
                              %%%          end;
        end;


        ACQEND = read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_ACQEND);
        if ( ACQEND ),


          % the last MDH and ADC line contains a 16-sample measurement of
          % unknown utility, so its stored in 'adc_cplx' without being
          % overwritten, but is not returned.

          DO__break_from_innerloop = 1;
          break;

        end;  % IF "ACQEND"


        %------------------------------------------------------------------------%
        % store measurement lines

        pos = [  1,                               % pos(01), columns
                 1 + mdh(idx).ushLine,            % pos(02)
                 1 + mdh(idx).ushChannelId,       % pos(03)
                 1 + mdh(idx).ushSet,             % pos(04)
                 1 + mdh(idx).ushEcho,            % pos(05)
                 1 + mdh(idx).ushPhase,           % pos(06)
                 1 + mdh(idx).ushRepetition,      % pos(07)
                 1 + mdh(idx).ushSeg,             % pos(08)
                 1 + mdh(idx).ushPartition,       % pos(09)
                 1 + mdh(idx).ushSlice,           % pos(10)
                 1 + mdh(idx).ushIda,             % pos(11)
                 1 + mdh(idx).ushIdb,             % pos(12)
                 1 + mdh(idx).ushIdc,             % pos(13)
                 1 + mdh(idx).ushIdd,             % pos(14)
                 1 + mdh(idx).ushIde,             % pos(15)
                 1 + mdh(idx).ushAcquisition      % pos(16)  % note: acquisition is same as average
              ];

        if ( DO__EXTRACT_REPETITION ),
          % shift repetition loopcounter
          pos(07) = pos(07) - iRep + 1;
        end;


        % (all this is to cater to the 7T host's peculiarities. sigh...)
        % this is also helpful for VD-platform data, where channel list is non-contiguous
        if ( ~FLAG__stored_channel_indices ),
          if ( ~ismember(pos(03), channel_indices) ),
            % accumulate list of channel indices
            channel_indices(end+1) = pos(03);

            % record / update first and last channel numbers
            channel_1 = channel_indices(1) - 1;
            channel_N = channel_indices(end) - 1;
          else,
            % all channels have been stored
            FLAG__stored_channel_indices = 1;

            % weirdness found
            if ( DO__SQUEEZE_CHANNELS && ...
                 ( (channel_indices(1) ~= 1) || ...
                   (channel_indices(end) ~= length(channel_indices)) )   ),
              FLAG__channel_squeeze = 1;

              % the coils should at LEAST be in order even if they are not contiguous
              [channel_sort, channel_pos] = sort(channel_indices);

              if ( ~isequal(channel_sort, channel_indices) ),
                warning('numerical receiver channel ID improperly ordered, unsorted!')
              end;

              channel_index_map = zeros(1,max(channel_indices));
              channel_index_map(channel_sort) = channel_pos;

              % after all the coil channels have been recorded, we need to go
              % back and resort the data that's been stored out-of-order.
              % since we don't know which data has been seen so far, test
              % everything. (this will only happen once.)
              if ( DO__RETURN_CELL_ARRAY ),
                if ( ~isequal(meas.data,{[]}) ), meas.data = meas.data(:,:,channel_sort,:,:,:,:,:,:,:,:,:,:,:,:,:); end;
              else,
                if ( ~isempty(meas.data) ), meas.data = meas.data(:,:,channel_sort,:,:,:,:,:,:,:,:,:,:,:,:,:); end;
              end;

              if ( ~isempty(timestamp_data) ),
                timestamp_data =           timestamp_data(          :,:,channel_sort,:,:,:,:,:,:,:,:,:,:,:,:,:);
              end;
              if ( ~isempty(timestamp_data_phascor1d) ),
                timestamp_data_phascor1d = timestamp_data_phascor1d(:,:,channel_sort,:,:,:,:,:,:,:,:,:,:,:,:,:);
              end;


              if ( FLAG__data_phascor1d ),     data_phascor1d     = data_phascor1d(       :,:,channel_sort,:,:,:,:,:,:,:,:,:,:,:,:,:); end;
              if ( FLAG__data_phascor2d ),     data_phascor2d     = data_phascor2d(       :,:,channel_sort,:,:,:,:,:,:,:,:,:,:,:,:,:); end;
              if ( FLAG__data_fieldmap ),      data_fieldmap      = data_fieldmap(        :,:,channel_sort,:,:,:,:,:,:,:,:,:,:,:,:,:); end;
              if ( FLAG__offline && DO__SAVE_OFFLINE ),  offline  = offline(              :,:,channel_sort,:,:,:,:,:,:,:,:,:,:,:,:,:); end;
              if ( FLAG__noiseadjscan ),       noiseadjscan       = noiseadjscan(         :,:,channel_sort,:,:,:,:,:,:,:,:,:,:,:,:,:); end;
              if ( FLAG__patrefscan ),         patrefscan         = patrefscan(           :,:,channel_sort,:,:,:,:,:,:,:,:,:,:,:,:,:); end;
              if ( FLAG__patrefscan_phascor ), patrefscan_phascor = patrefscan_phascor(   :,:,channel_sort,:,:,:,:,:,:,:,:,:,:,:,:,:); end;

              % hack due to logic used with FLAG__smsrefscan
              if ( ~isempty(smsrefscan) ),     smsrefscan         = smsrefscan(           :,:,channel_sort,:,:,:,:,:,:,:,:,:,:,:,:,:); end;

              if ( FLAG__smsrefscan_phascor ), smsrefscan_phascor = smsrefscan_phascor(   :,:,channel_sort,:,:,:,:,:,:,:,:,:,:,:,:,:); end;
              if ( FLAG__refphasestabscan ),   refphasestabscan   = refphasestabscan(     :,:,channel_sort,:,:,:,:,:,:,:,:,:,:,:,:,:); end;
              if ( FLAG__phasestabscan ),         phasestabscan   =    phasestabscan(     :,:,channel_sort,:,:,:,:,:,:,:,:,:,:,:,:,:); end;
              if ( FLAG__phasestabscan ),   patrefphasestabscan   = patrefphasestabscan(  :,:,channel_sort,:,:,:,:,:,:,:,:,:,:,:,:,:); end;

            end;
          end;
        end;


        if ( DO__CANONICAL_REORDER_COIL_CHANNELS ),
          if ( channel_1 ~= 0 ),
            disp_optional('non-contiguous channels found---aborting canonical reordering.');
            DO__CANONICAL_REORDER_COIL_CHANNELS = 0;
          else,
            % cannot reorder canonically and squeeze channels, so if by this
            % point the canonical reordering is still valid, disable squeeze
            % channels
            FLAG__channel_squeeze = 0;
          end;
          pos(03)  = coil_order(pos(03));
        end;

        % apply map
        if ( FLAG__channel_squeeze ),
          %if ( ismember(pos(03), channel_indices) ),
          if ( any(pos(03) == channel_indices) ),
            pos(03) = channel_index_map( pos(03) );
          end;
        end;

        % if an SMS scan has already been detected
        if ( DO__SMSREFSCAN && ~FLAG__smsrefscan_completed && FLAG__smsrefscan ),

          % if this slice has been visited, then SMS reference lines have been stored.
          % (NOTE: this logic assumes that SMS reference lines appear before
          % all other data categories that use the slice loopcounter)
          if ( [sms_slc_buffer >= 0 ] && [pos(10) < sms_slc_buffer] && [pos(10) == 1] ),
            FLAG__smsrefscan_completed = 1;
          end;

        end;


        % correct for idiosyncratic, non-contiguous indexing schemes employed by
        % Siemens [found in all versions (starting from VB11A) examined]

        % Segment: all lines acquired in the negative readout direction (labeled
        % with the MDH_REFLECT flag) are assigned to their own segment via the
        % loopcounter to be compatible with the "OnlineTSE" phase correction
        % functor. as a result, a multishot segmented acquisition with
        % 'NSegMeas' legitimate segments will have 2*NSegMeas segments in the
        % loop counter---the even segments will be normal lines and the odd
        % segments will be reflected lines (0-based indexing). we must then
        % collapse the even and odd lines from each shot into a single segment.
        % [TODO: fix this hack...not all sequences use this scheme]

        % first, check if this line is a PHASCOR line (in case we're at the first line)
        if ( ~FLAG__phascor1d && read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_PHASCOR) ),
          FLAG__phascor1d = 1;
        end;

        % ASSUME that if any phascor lines have been collected that all lines
        % (e.g., data and "patrefscan" lines) follow the same convention with
        % the segment loop counter---but this may not be true!
        if ( FLAG__phascor1d && FLAG__data_reflect && DO__PHASCOR_COLLAPSE_SEGMENTS ),
          pos(08) = pos(08)/2;  % if condition true, pos(08) should always be even, so quotient should be integer!
        end;

        % Acquisition: the last of reference navigators scan batch [by
        % default three are collected for each slice (2D) or partition (3D)]
        % is assigned an acquisition index of 1, for siemens's
        % "EPIPhaseCorrPEFunctor" [TODO: fix this hack]
        if ( read_meas__extract_bit(aulEvalInfoMask, ...
                                    EvalInfoMask.MDH_PHASCOR) ),
          pos(16) = ceil(pos(16)/2);
        end;


        %%% begin logic to determine category of each ADC line

        %%% CATEGORY #1: noise adjustment scan lines
        if ( read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_NOISEADJSCAN) ),

          if ( ~FLAG__noiseadjscan ),
            disp_optional(' NOISEADJSCAN detected.');
            FLAG__noiseadjscan = 1;
            noise_line = 0;
          end;

          % increment noise adjust scan line counter after all channels' data is in
          if ( mdh(idx).ushChannelId == channel_1 ),
            noise_line = noise_line + 1;
          end;

          % TODO: establish line counter for noise scans

          % remove fft scaling of noise
          if ( DO__APPLY_FFT_SCALEFACTORS ),
            adc_cplx = adc_cplx / fft_scale( double(mdh(idx).ushChannelId+1) );
          end;
          noiseadjscan(:, noise_line, pos(03)) = adc_cplx;

          continue;

          %%% CATEGORY #3: phase stabilization scans
        elseif ( read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_PHASESTABSCAN) ),

          % HACK: for now, any line that is flagged as a PHASESTABSCAN,
          % regardless if it is an image line or a ACS line, will be
          % accumulated in the same array. this is because, so far, i have
          % only encountered both phase stabilization and iPAT in scans where
          % the ACS lines are integrated (i.e., the ACS lines are
          % PATREFANDIMASCAN lines). if the ACS lines are acquired separately,
          % as they are in EPI time series data for example, there would have
          % to be a separate array to accumulate the phase stabilization
          % navigators for the ACS lines in addition to the array that
          % accumulates the phase stabilization navigators for the image lines
          % (analogous to the phase correction navigators and the
          % "patrefscan_phascor" array).

          if ( ~FLAG__phasestabscan && ...
               ~read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_REFPHASESTABSCAN)),

            % need the time of the phase stabilization navigator echo, "TS",
            % for image reconstruction, but its NOT available from the header
            % for some weird reason. extract it from the MDH header. it should
            % be the same across segments and slices (which are the only
            % dimensions over which its repeated).
            TS = double(mdh(idx).ulTimeSinceLastRF) / 1e3;    % store as milliseconds

            if ( TS == TS0 ),
              disp_optional(        ' PHASESTABSCAN detected.');
            else,
              disp_optional(sprintf(' PHASESTABSCAN detected;     (( TS = %.3f ms )).', TS));
            end;
            FLAG__phasestabscan = 1;

          end;

          % NOTE: both the phase stabilization measurements and the phase
          % stabilization reference measurements are stored as echo 0
          % regardless of how many echoes are acquired for the imaging data.

          %%% VB13A ICE manual, V0.7, p. 274:
          %  "Phase stabilization is restricted to the sharing of phase
          %   stabilization scans between different contrasts (echoes), i.e.
          %   the imaging scans, which belong to the same phase stabilization
          %   echo, may only differ in the ECO ICE Dimension."


          %%% CATEGORY #3A: phase stabilization reference
          if ( read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_REFPHASESTABSCAN) ),


            if ( ~FLAG__refphasestabscan ),

              % need the time of the phase stabilization navigator echo, "TS",
              % for image reconstruction, but its NOT available from the header
              % for some weird reason. extract it from the MDH header. it should
              % be the same across segments and slices (which are the only
              % dimensions over which its repeated).
              TS0 = double(mdh(idx).ulTimeSinceLastRF) / 1e3;    % store as milliseconds

              disp_optional(sprintf(' REFPHASESTABSCAN detected;  (( TS = %.3f ms )).', TS0));
              FLAG__refphasestabscan = 1;

            end;

            % the phase stabilization works by collecting a full set
            % (potentially multiple echoes) of lines at the beginning of each
            % slice that are tagged as MDH_REFPHASESTABSCAN, followed by a
            % *single* measurement that is tagged as both MDH_REFPHASESTABSCAN
            % and MDH_PHASESTABSCAN. it is this last scan that serves as the
            % reference, and it will exhibit the same timing as the subsequent
            % phase stabilization scans that follow each group of image
            % echoes. therefore the group of MDH_REFPHASESTABSCAN lines
            % preceding each *single* reference echo are not used by the
            % correction and are discarded.

            refphasestabscan(...
                :, pos(02), ...
                pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) ...
                = adc_cplx;

            refphasestabtime(...
                1, pos(02), ...
                pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) ...
                = mdh(idx).flTimeStamp_ms;
          else,

            %%% CATEGORY #3B: phase stabilization navigators for ACS data
            if ( read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_PATREFANDIMASCAN) ),

              patrefphasestabscan(...
                  :, pos(02), ...
                  pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                  pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) ...
                  = adc_cplx;

              patrefphasestabtime(...
                  1, pos(02), ...
                  pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                  pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) ...
                  = mdh(idx).flTimeStamp_ms;

              phasestabscan(...
                  :, pos(02), ...
                  pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                  pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) ...
                  = adc_cplx;

              phasestabtime(...
                  1, pos(02), ...
                  pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                  pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) ...
                  = mdh(idx).flTimeStamp_ms;

              %%% CATEGORY #3C: phase stabilization navigators for ACS & image data
            elseif ( read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_PATREFSCAN) ),

              patrefphasestabscan(...
                  :, pos(02), ...
                  pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                  pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) ...
                  = adc_cplx;

              patrefphasestabtime(...
                  1, pos(02), ...
                  pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                  pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) ...
                  = mdh(idx).flTimeStamp_ms;

              %%% CATEGORY #3D: phase stabilization navigators for image data
            else,

              phasestabscan(...
                  :, pos(02), ...
                  pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                  pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) ...
                  = adc_cplx;

              phasestabtime(...
                  1, pos(02), ...
                  pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                  pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) ...
                  = mdh(idx).flTimeStamp_ms;

            end;


          end;
          continue;

          %%% CATEGORY #2: iPAT ACS lines reference scan
        elseif ( read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_PATREFSCAN) ),

          if ( ~FLAG__patrefscan ),
            disp_optional(' PATREFSCAN detected.');
            FLAG__patrefscan = 1;

            pat_acq_counter     = pos(16);

%            %% --> BEGIN geste mod, 2013/aug/26
%            pat_acq_counter_new = pos(16);
%            pat_acq_counter_sign = -1;
%            %% <-- END geste mod, 2013/aug/26

          end;

          %% --> BEGIN geste mod, 2013/oct/26
          if ( DO__TENC_ACS ),
	    pos(16) = double(mdh(idx).aushFreePara(1)) + 1;
	  end;
          %% <-- END geste mod, 2013/oct/26


          % works with PHASCOR lines assuming that pos(16) is divided by 2!
          if ( pos(16) > pat_acq_counter ),
            if ( ~FLAG__repeated_patrefscan ),
              % detected Benner's implementation of multiple ACS scan averages
              disp_optional('  multiple PATREFSCAN acquistions detected.');
              FLAG__repeated_patrefscan = 1;
            end;

            disp_optional(sprintf('   acq: %02d', pat_acq_counter));
            pat_acq_counter = pos(16);
          end;


          %%% CATEGORY #2A: phase correction navigator line for iPAT ACS lines reference scan
          if ( read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_PHASCOR) ),

            if ( ~FLAG__patrefscan_phascor ),
              disp_optional(' PATREFSCAN_PHASCOR detected.');
              FLAG__patrefscan_phascor = 1;

              % check if first phase correction navigator line is reflected (see comments below)
              if ( read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_REFLECT) ),
                FLAG__patrefscan_phascor_orderswap = 1;
              end;

              % this shouldn't be necessary...but the MDH_FIRSTSCANINSLICE
              % flag doesn't get set for ACS lines! grr... pat_nav_slice is
              % used as memory: it tracks the last pat_nav slice. if current
              % slice number is greater than last slice number, a new slice
              % has been encountered. similarly, if the current slice is
              % less than the last slice number, a new cycle through the
              % slices is encountered

              pat_nav_slice = 0;
              pat_nav_base_segment = pos(08);
              pat_nav_line = 0;
              pat_nav_rep = pos(07);
              pat_nav_seg = 0;

            end;

            % reset navigator scan line counter at the beginning of each slice
            if (  (mdh(idx).ushChannelId == channel_1) && ...
                  ( read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_FIRSTSCANINSLICE) || ...
                    pos(10) > pat_nav_slice )  ),
              pat_nav_line = 0;
            end;

            % if this slice has been seen before, reset line counter
            if ( pos(10) < pat_nav_slice ),
              pat_nav_line = 0;

              % if it is because it is a new repetition, also reset segment counter
              if ( pos(07) > pat_nav_rep ),
                pat_nav_rep = pos(07);
                pat_nav_seg = 0;
              else,
%               disp(sprintf('%2d', pos));
%                pat_nav_seg = pat_nav_seg+2;
              end;

              % NOTE: not sure how to differentiate between new shots and
              % multiple repetitions...

            end;

            % reset navigator scan line counter (rare---occurs only for TSE data?)
            if ( (mdh(idx).ushChannelId == channel_1) && ...
                 pos(08) >= (pat_nav_base_segment+2) )
              pat_nav_line = 0;
              pat_nav_base_segment = pos(08);
              disp('pat_nav_line reset #3');
            end;

            % increment navigator scan line counter after all channels' data is in
            if ( mdh(idx).ushChannelId == channel_1 ),
              pat_nav_line = pat_nav_line + 1;
            end;

            dims = size(patrefscan_phascor);
            dims(end+1:16) = 1;

            %% --> BEGIN geste mod, 2013/aug/26

%            pat_acq_thisscan = [...
%                  1,   pat_nav_line, ...
%                  pos(03), pos(04), pos(05), pos(06), pat_nav_rep, pat_nav_seg+pos(08), 1, ...
%                  pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)];
%
%            pat_acq_thisscan_mod = pat_acq_thisscan;
%            if ( pat_acq_counter_sign < 0 ),
%              pat_acq_thisscan_mod(08) = pat_acq_thisscan_mod(08) + 1;
%            end;
%
%            if ( isequal(pat_acq_thisscan_mod, pat_acq_firstscan) ),
%              scan_num
%              pat_acq_thisscan
%              pat_acq_counter_sign = pat_acq_counter_sign * (-1);
%              pat_acq_counter_new = pat_acq_counter_new + pat_acq_counter_sign;
%            end;
%
%            pos(16) = pat_acq_counter_new;
%
            %% <-- END geste mod, 2013/aug/26



            % attempt to store an array of visited k-space locations to
            % check whether any incoming lines have loopcounters that would
            % lead to clobbering/overwritting previously stored data---which
            % could be caused either by bad logic here or by a bug in the
            % sequence

%           visited = 0;
%            try,
%             if ( any(...
%                 patrefscan_phascor(...
%                     :,   pat_nav_line, ...     % only center line and center partition collected
%                     pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), 1, ...
%                     pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16))) ),
%               visited = 1;
%             end;
%           catch,
              patrefscan_phascor(...
                  :,   pat_nav_line, ...     % only center line and center partition collected
                  pos(03), pos(04), pos(05), pos(06), pat_nav_rep, pat_nav_seg+pos(08), 1, ...
                  pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = adc_cplx;
%           end;

           if ( isempty(pat_acq_firstscan) ),
             pat_acq_firstscan = [...
                  1,   pat_nav_line, ...
                  pos(03), pos(04), pos(05), pos(06), pat_nav_rep, pat_nav_seg+pos(08), 1, ...
                  pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)];
           end;


%           if ( visited ), warning(sprintf('visited line detected: scan number %d', scan_num)); end;

            % this STILL shouldn't be necessary...but the MDH_FIRSTSCANINSLICE flag doesn't get set for ACS lines! grr...
            pat_nav_slice = pos(10);


            %%% CATEGORY #2B: lines used for both iPAT ACS and data
          elseif ( read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_PATREFANDIMASCAN) ),
            % there are three PAT reconstruction modes: Integrated, Separate,
            % and Averaging. in 'Integrated' mode, the ACS lines are acquired
            % in-line with the data, so some serve both as data and ACS, whereas
            % in 'Separate' mode the ACS lines are acquired separately (e.g., as
            % in accelerated EPI). ('Averaging' mode is not supported
            % here....yet.)

            % for 'Integrated' mode, the ACS and data lines share the same
            % loopcounter numbering scheme, so the direct approach to storing
            % the k-space data would yield an inefficient storage where omitted
            % lines of k-space are stored as zeros! for now, let's be lazy and
            % waste some memory since this mode is probably just for
            % single-acquisition structural (e.g., 3D-MPRAGE) images, so we
            % don't have to be clever in storing the normal data lines.

            % NOTE: in this case, the last skipped k-space lines will not be
            % present, so the iPAT data will contain R-1 fewer lines than the
            % corresponding non-accelerated data. we could just append those R-1
            % lines of zeros to the end, but that would require knowing the
            % value of R. hmmm... nah, skip it.

            % so if a line is both ACS and an image line, store it *twice*: once
            % as data and once as reference.

            % (for 1-dimensional acceleration, R could be estimated from the
            % data sparsity calculated at the end, but this would not work for
            % 2-dimensional imaging!)

            if ( ~FLAG__patrefandimascan ),
              disp_optional(' PATREFANDIMASCAN detected.');
              FLAG__patrefandimascan = 1;
              FLAG__data = 1;
            end;

            if ( DO__RETURN_CELL_ARRAY ),   % 1 of 4
              patrefscan{...
                :, pos(02), ...
                pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)} = adc_cplx;
            else,
              patrefscan(...
                :, pos(02), ...
                pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = adc_cplx;
            end;

            if ( pos(02) > maxLin ),
              maxLin = pos(02);
              if ( DO__DISP_PROGRESS == 2 ), fprintf('  <%s>  l: %04d\n', datestr(now, 13), maxLin); end;
            end;


            % store ADC line in listed position within data volume (using 1-based indices)
            if ( ~DO__RETURN_PATREF_ONLY ),
              if ( DO__RETURN_CELL_ARRAY ),   % 2 of 4
                meas.data{...
                    :, pos(02), ...
                    pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                    pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)} ...
                    = adc_cplx;
              else,
                % hack, hack, hack
                if ( isempty(meas.data) ),
                  meas.data = complex(zeros(size(adc_cplx), matrixtypestr));
                end;

                meas.data(...
                    :, pos(02), ...
                    pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                    pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) ...
                    = adc_cplx;
              end;
            end;

            %%% CATEGORY #2C: just an ordinary ACS line
          else,

%            %% --> BEGIN geste mod, 2013/aug/26
%            pos(16) = pat_acq_counter_new;
%            %% <-- END geste mod, 2013/aug/26

            if ( DO__RETURN_CELL_ARRAY ),   % 3 of 4
              patrefscan{...
                :, pos(02), ...
                pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)} = adc_cplx;
            else,
              patrefscan(...
                :, pos(02), ...
                pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = adc_cplx;
            end;

            %          scan_patrefscan(...
            %              1, pos(02), ...
            %              pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
            %              pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = mdh(idx).ulScanCounter;

            if ( pos(02) > maxLin ),
              maxLin = pos(02);
              if ( DO__DISP_PROGRESS == 2 ), fprintf('  <%s>  l: %04d\n', datestr(now, 13), maxLin); end;
            end;

          end;
          continue;

          %%% CATEGORY #5: SMS reference data
        elseif ( DO__SMSREFSCAN && ~FLAG__smsrefscan_completed ),

          if ( ~FLAG__smsrefscan ),
            % the SMS reference isn't actually detected, but is defined as the
            % first repetition!
            disp_optional(' SMSREFSCAN mode *selected by user*.');
            FLAG__smsrefscan = 1;

            sms_acq_counter = pos(16);

            % count slices, and when all slices are in the SMS reference data is completed
            % (NOTE: this scheme doesn't allow for multiple acquisitions!)
            sms_slc_buffer = -1;

          end;

          % works with PHASCOR lines assuming that pos(16) is divided by 2!
          if ( pos(16) > sms_acq_counter ),
            if ( ~FLAG__repeated_smsrefscan ),
              % for a future implementation of multiple ACS scan averages
              disp_optional('  multiple SMSREFSCAN acquistions detected.');
              FLAG__repeated_smsrefscan = 1;
            end;

            disp_optional(sprintf('   acq: %02d', sms_acq_counter));
            sms_acq_counter = pos(16);
          end;

          %%% CATEGORY #5A: phase correction navigator line for SMS reference scan
          if ( read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_PHASCOR) ),

            if ( ~FLAG__smsrefscan_phascor ),
              disp_optional(' SMSREFSCAN_PHASCOR detected.');
              FLAG__smsrefscan_phascor = 1;

              % check if first phase correction navigator line is reflected (see comments below)
              if ( read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_REFLECT) ),
                FLAG__smsrefscan_phascor_orderswap = 1;
              end;

              % perhaps the SMS reference lines will in the future behave like ACS reference lines?
              sms_nav_slice = 0;

              sms_nav_base_segment = pos(08);

              sms_nav_line = 0;

            end;

            % reset navigator scan line counter at the beginning of each slice
            if (  (mdh(idx).ushChannelId == channel_1) && ...
                  ( read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_FIRSTSCANINSLICE) || ...
                    pos(10) > sms_nav_slice )  ),
              sms_nav_line = 0;
            end;

            if ( pos(10) < sms_nav_slice ),
              sms_nav_line = 0;
            end;


            % reset navigator scan line counter
            if ( (mdh(idx).ushChannelId == channel_1) && ...
                 pos(08) >= (sms_nav_base_segment+2) )
              sms_nav_line = 0;
              sms_nav_base_segment = pos(08);
            end;


            % increment navigator scan line counter after all channels' data is in
            if ( mdh(idx).ushChannelId == channel_1 ),
              sms_nav_line = sms_nav_line + 1;
            end;

            smsrefscan_phascor(...
                :,   sms_nav_line, ...     % only center line and center partition collected
                pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), 1, ...
                pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = adc_cplx;


            % again, using the same logic for iPAT lines for the SMS lines
            sms_nav_slice = pos(10);
          else,

            smsrefscan(...
                :, pos(02), ...
                pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = adc_cplx;

            % record slices that are stored as SMS reference scans
            if ( pos(10) > sms_slc_buffer ),
              sms_slc_buffer = pos(10);
            end;

          end;

          continue;

          %%% CATEGORY #4: data
        else,

          if ( DO__EXTRACT_REPETITION && (firstdata_pos < 1) ),
            firstdata_pos = meas_pos;
            fseek(fp, meas_pos, 'bof');
            continue;
          end;


          %%% CATEGORY #4A: phase correction navigator line for data
          if ( read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_PHASCOR) ),

            if ( ~FLAG__data_phascor1d ),
              disp_optional(' PHASCOR (1D) detected.');
              FLAG__data_phascor1d = 1;

              if ( data_begin < 1 ),
                data_begin = meas_pos;
              end;

              new_seg = -1;
              new_slc = -1;
              % check if first phase correction navigator line is reflected
              if ( read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_REFLECT) ),

                % for example, siemens's "ep2d_bold" reflects the first
                % navigator line, in which case the EVEN-indexed navigator
                % lines will correspond to the ODD-indexed data lines. yuck!
                FLAG__data_phascor1d_orderswap = 1;

                % TODO: store the reflected and non-reflected lines in
                % separate matrices, instead of adhering to the "odd--even"
                % organization, since siemens appears to like to collect more
                % of one than the other! (DONE!)
              end;

              nav_base_segment = pos(08);

              segment = 0;

              % initialize
              nav_line = 0;

            end;


            % reset navigator scan line counter at the beginning of each slice
            if ( mdh(idx).ushChannelId == channel_1 && ...
                 read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_FIRSTSCANINSLICE) ),
              nav_line = 0;
              nav_reset = nav_line;
              segment = 1;
              nav_base_segment = pos(08);
            end;

            %%%       old_seg = new_seg;
            %%%       new_seg = mdh(idx).aushFreePara(2) + 1;
            %%%
            %%%       old_slc = new_slc;
            %%%       new_slc = pos(10);
            %%%
            %%%       if ( mdh(idx).ushChannelId == channel_1 && (new_seg ~= old_seg) ),
            %%%         disp('new segment!');
            %%%         segment = segment + 1;
            %%%
            %%%         % note: this is redundant, since new segments and new slices will both occur at same MDHs
            %%%   %     nav_line = 0;
            %%%       end;
            %%%
            %%%       % because later segments do not use FIRSTSCANINSLICE flag, we need
            %%%             % some memory and detect whether the slice has changed
            %%%       if ( mdh(idx).ushChannelId == channel_1 && (new_slc ~= old_slc) ),
            %%%   %     segment = 1;
            %%%               nav_line = 0;
            %%%       end;
            %%%
            %%%       pos(04) = segment;


            % if a truly a new segment is found, reset navigator scan line counter
            % (happens with TSE)
            if ( (mdh(idx).ushChannelId == channel_1) && ...
                 pos(08) >= (nav_base_segment+2) )
              nav_line = nav_reset;
              nav_base_segment = pos(08);
            end;

            %         pos(08) = pos(08) + segment;
            %
            %         pos_virtual = double(pos.');
            %          pos_virtual(02) = max([nav_line,1]);
            %          pos_virtual(09) = 1;
            %
            %          % HACK for conventional (slice-interleaved) multi-shot
            %          if ( ~isempty(data_phascor1d) ),
            %            fail = 0;
            %           stored = 0;
            %            try,
            %              stored = data_phascor1d(pos_virtual(01), pos_virtual(02), pos_virtual(03), ...
            %                             pos_virtual(04), pos_virtual(05), pos_virtual(06), ...
            %                             pos_virtual(07), pos_virtual(08), pos_virtual(09), ...
            %                             pos_virtual(10), pos_virtual(11), pos_virtual(12), ...
            %                             pos_virtual(13), pos_virtual(14), pos_virtual(15), ...
            %                             pos_virtual(16));
            %            catch,
            %              fail = 1;
            %            end;
            %
            %            if ( ~fail && stored ~= 0 ),
            %
            %             nav_line = 0;
            %             disp('--------------------------------------------');
            %             segment = segment + 2
            %             pos(08) = pos(08) + 2;
            %             disp('--------------------------------------------');
            %
            %            end;
            %          end;


            % increment navigator scan line counter after all channels' data is in
            if ( mdh(idx).ushChannelId == channel_1 ),
              nav_line = nav_line + 1;
            end;

            if ( DO__EXTRACT_REPETITION && ~FLAG__FOUND_EXTRACT_REPETITION )
              continue;
            end;

            % if 2D, only one partition. if 3D each navigator line is
            % collected at center of k-space, but partition loopcounter set
            % to center partition, so here we force it to 1 so that only
            % nav_line counter is incremented.

            data_phascor1d(...
                :, nav_line, ...     % only center line and center partition collected
                pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), 1, ...
                pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = adc_cplx;

            if ( DO__STORE_TIMESTAMP ),
              timestamp_data_phascor1d(...
                  :, nav_line, ...     % only center line and center partition collected
                  pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), 1, ...
                  pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = mdh(idx).flTimeStamp_ms / 1e3;
            end;

            % navigator scan partition counter should NEVER be reset


            % PE line and PE partition incrementing depends on whether
            % acquisition is 2D or 3D. so can check whether line and partition
            % loop counter values on first line are equal to the center
            % line/partition of k-space.

            %          % increment navigator scan partition counter at the end of each partition
            %          if ( FLAG__stored_channel_indices && ...
            %               ( mdh(idx).ushChannelId == channel_N ) && ...
            %               read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_LASTSCANINSLICE) ...
            %               && ( mdh(idx).ushPartition == mdh(idx).ushKSpaceCentrePartitionNo ) ...
            %               && ( mdh(idx).ushKSpaceCentrePartitionNo > 0 ) ...    % since always == 0 for 2D imaging
            %               ),
            %            nav_part = nav_part + 1;
            %          end;

            % QUESTION: for two partitions (pathological case), is the center == 0 or == 1?


            %%%/        if (  ( (mdh(idx).ushChannelId + 1) == mdh(idx).ushUsedChannels ) && ...
            %%%/              ( (mdh(idx).ushAcquisition + 1) > 1 )  ),
            %%%/
            %%%/          % for some multi-shot scans, get navigator
            %%%/          % lines with each shot
            %%%/          nav_ = nav_ + 1;
            %%%/        end;

            %%% CATEGORY #4B: field map lines for data
            % (PHASCOR2D HACK: currently the only signifier for 2D phascor data is
            % the absence of the MDH_ONLINE and MDH_PHASCOR flags)
          elseif ( ~read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_ONLINE) ),

            if ( read_meas__extract_bit(aulEvalInfoMask, EvalInfoMask.MDH_PHASCOR) ),

              if ( ~FLAG__data_phascor2d ),
                disp_optional(' PHASCOR (2D) detected.');
                FLAG__data_phascor2d = 1;
                nav_part = 1;
              end;

              %%%%      if ( mdh(idx).ushChannelId == channel_N ),
              %%%%        fprintf('     [%s]  %04d\n', datestr(now, 13), scan_num);
              %%%%      end;

              % navigator scan partition counter should NEVER be reset


              %% [2011/sep/16] do not need these pesky OFFLINE lines any
              %% more; in future will % add an option to allow user to request
              %% the OFFLINE scans which are often used when prototyping
              %% sequences that don't yet have working ICE code.

              data_phascor2d(...
                  :, pos(02), ...
                  pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), nav_part, ...  % only center partition collected
                  pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = adc_cplx;

              %         scan_phascor2d(...
              %              1, pos(02), ...
              %              pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), nav_part, ...
              %              pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = mdh(idx).ulScanCounter;


              % increment navigator scan partition counter at the end of each partition (signified with the "LASTSCANINSLICE" flag)
              % (ASSUME here that these navigator lines ONLY appear in 3D sequences)
              if ( FLAG__stored_channel_indices && ...
                   ( mdh(idx).ushChannelId == channel_N ) ...
                   && read_meas__extract_bit(aulEvalInfoMask, ...
                                             EvalInfoMask.MDH_LASTSCANINSLICE) ),
                nav_part = nav_part + 1;
              end;

              %%% CATEGORY #4C: "offline" lines
            else,

              if ( ~FLAG__offline ),
                disp_optional(' OFFLINE detected.');
                FLAG__offline = 1;
              end;

              if ( DO__SAVE_OFFLINE ),
                offline(...
                    :, pos(02), ...
                    pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                    pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) = adc_cplx;
              end;
            end;

            %%% CATEGORY #4D: data lines (finally!)
          else,

            if ( ~FLAG__data ),
              FLAG__data = 1;
            end;

	    if ( DO__RETURN_PHASCOR_ONLY ),
	      break;              % break from FOR loop
	    end;
            
	    
	    
	    % store ADC line in listed position within data volume (using 1-based indices)
            if ( ~DO__RETURN_PATREF_ONLY ),
              if ( DO__RETURN_CELL_ARRAY ),   % 4 of 4
                meas.data{...
                    :, pos(02), ...
                    pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                    pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)} ...
                    = adc_cplx;
              else,

                % hack, hack, hack
                if ( isempty(meas.data) ),
                  meas.data = complex(zeros(size(adc_cplx), matrixtypestr));
                end;

		%                                                    2    4    5    6    7    8    9   0   
		if ( DO__DISP_PROGRESS > 2 &&  mdh(idx).ushChannelId == channel_N ), 
		  fprintf('  <%s>  n: %04d %04d %04d %04d %04d %04d %04d %04d\n', datestr(now, 13), pos(02), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), pos(10));
		end;
		
		% TODO: SelectChannels should exclude SqueezeChannels and CanonicalReorderCoilChannels
		if ( DO__SELECT_CHANNELS && ~ismember(pos(03), coil_selected) ),
		  continue;
		end;
		
		if ( DO__SELECT_CHANNELS ),
		  pos(03) = find(pos(03) == coil_selected);
		end;
		
		meas.data(...
                    :, pos(02), ...
                    pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                    pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) ...
                    = adc_cplx;
              end;
            else,
              if ( FLAG__patrefscan && pos(07) >= 2 ),
                % assume that repetitions are only used with iPAT for
                % acquisitions (like EPI) where reference scans are
                % separate
                DO__break_from_innerloop = 1;
                break;
              end;
            end;
            if ( DO__STORE_TIMESTAMP ),
              timestamp_data(...
                  :, pos(02), ...
                  pos(03), pos(04), pos(05), pos(06), pos(07), pos(08), pos(09), ...
                  pos(10), pos(11), pos(12), pos(13), pos(14), pos(15), pos(16)) ...
                  = mdh(idx).flTimeStamp_ms / 1e3;
            end;


            if ( pos(07) > maxRep ),
              maxRep = pos(07);
              if ( DO__DISP_PROGRESS ), fprintf('  <%s>  r: %04d\n', datestr(now, 13), maxRep); end;
            end;

            if ( pos(02) > maxLin ),
              maxLin = pos(02);
              if ( DO__DISP_PROGRESS == 2 ), fprintf('  <%s>  l: %04d\n', datestr(now, 13), maxLin); end;
            end;


            if ( mdh(idx).ushChannelId == channel_N ),
              %              fprintf('     [%s]  %04d\n', datestr(now, 13), scan_num);
            end;


          end;
        end;  % IF (to determine category assignment)
      end;  % FOR (loop over channels in VD-platform data)
      %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
    end;  % WHILE (main loop)
    %======================================================================%
    % M A I N    L O O P        E N D
    %======================================================================%

  catch,

    caught = lasterror;
    status = dbstatus;

    disp(sprintf('<!> [%s]:  trouble at line %d...', mfilename, caught.stack(1).line));
%    disp(lasterr);


    if ( feof(fp) ),

      caught.message = sprintf('EOF encountered before ACQEND! last scan number = %d ', ...
                               scan_num);

      if ( DO__RECOVER_FROM_INCOMPLETE ),

        fprintf(1, '\n');
        warning(caught.identifier, caught.message);
        fprintf(1, '\n');

        disp(sprintf('<!> [%s]:  ABORTING read, incomplete file!!!', mfilename));
        fprintf(1, '\n');

        FLAG__status_OK = 0;

      else,
        error(caught.identifier, caught.message);
      end;

    else,

      if ( 0 ),
        rethrow(caught);
      else,
        err = lasterr;
        disp(['??? ', err]); beep; pause(0.1); beep;
        keyboard
      end;


    end;


  end;


  %------------------------------------------------------------------------%
  %%% post file I/O reorganization

  % if the phase correction lines and the data lines mismatch in terms
  % whether the odd or even lines are reflected, shift the phase
  % correction line index by one and insert a duplicate of the second
  % line at the beginning.


  if ( DO__PHASCOR_COLLAPSE_SEGMENTS && FLAG__patrefscan_phascor_orderswap ),

    disp_optional('first "PATREFSCAN_PHASCOR" line reflected! compensating with odd--even order swap...');

    patrefscan_phascor(:, 2:end+1, :, :, :, :, :, :, :, :, :, :, :, :, :, :) = ...
        patrefscan_phascor(:, 1:end, :, :, :, :, :, :, :, :, :, :, :, :, :, :);

    patrefscan_phascor(:, 1, :, :, :, :, :, :, :, :, :, :, :, :, :, :) = ...
        patrefscan_phascor(:, 3, :, :, :, :, :, :, :, :, :, :, :, :, :, :);
  end;


  if ( DO__PHASCOR_COLLAPSE_SEGMENTS && FLAG__smsrefscan_phascor_orderswap ),

    disp_optional('first "SMSREFSCAN_PHASCOR" line reflected! compensating with odd--even order swap...');

    smsrefscan_phascor(:, 2:end+1, :, :, :, :, :, :, :, :, :, :, :, :, :, :) = ...
        smsrefscan_phascor(:, 1:end, :, :, :, :, :, :, :, :, :, :, :, :, :, :);

    smsrefscan_phascor(:, 1, :, :, :, :, :, :, :, :, :, :, :, :, :, :) = ...
        smsrefscan_phascor(:, 3, :, :, :, :, :, :, :, :, :, :, :, :, :, :);
  end;


  if ( DO__PHASCOR_COLLAPSE_SEGMENTS && FLAG__data_phascor1d_orderswap ),

    disp_optional('first "PHASCOR" line reflected! compensating with odd--even order swap...');

    data_phascor1d(:, 2:end+1, :, :, :, :, :, :, :, :, :, :, :, :, :, :) = ...
        data_phascor1d(:, 1:end, :, :, :, :, :, :, :, :, :, :, :, :, :, :);

    data_phascor1d(:, 1, :, :, :, :, :, :, :, :, :, :, :, :, :, :) = ...
        data_phascor1d(:, 3, :, :, :, :, :, :, :, :, :, :, :,: , :, :);
  end;


  %------------------------------------------------------------------------%

  fclose(fp);


  t1 = clock;
  runtime_seconds = etime(t1,t0);

  TIME = sprintf('[[ %02dh %02dm %02ds ]]', ...
                 fix(runtime_seconds/60/60), ...
                 rem(fix(runtime_seconds/60), 60), ...
                 rem(fix(runtime_seconds), 60));

  dstr = sprintf('total read time = %s', TIME);
  disp_optional(sprintf('<t> [%s]: %s', mfilename, dstr));

%
%  if ( nnz(data) == 0 ),
%    if ( nargout > 0 ),
%      warning('no data found!!!');
%    else,
%      disp_optional('no data found!!!');
%      return;
%    end;
%  end;
%
%
%  disp_optional(sprintf('data memory allocated:    %10.2f MB', getfield(whos('data'), 'bytes')/2^20));
%
%  density = nnz(data) / numel(data);
%  % sparsity = 1 - density;
%  redundancy = 1 / density;
%
%  if ( redundancy > 1 && density < 0.99 ),
%    disp_optional(sprintf('data sparsity detected---potential redundancy factor = [[ %5.1f X ]]', ...
%                 redundancy));
%  end;
%

  %------------------------------------------------------------------------%

  if ( DO__RETURN_STRUCT ),

    %meas = struct;
    meas.file = strcat(filestr, extstr);
    meas.readtime = datestr(t1, 'yyyy-mmm-dd HH:MM:SS');

    meas.acqtime = [meas_timestart, mdh(end).flTimeStamp_ms];

    % some functions (e.g., cell2mat) require that each entry has the same
    % datatype
    if ( DO__RETURN_CELL_ARRAY ),
      tic;
      disp(sprintf('==> [%s]: converting cell array to type "%s"', mfilename, matrixtypestr));
      meas.data = cellfun(matrixtype, meas.data, 'UniformOutput', 0);
      toc;
    end;

    if ( FLAG__data_phascor1d ),
      meas.data_phascor1d = data_phascor1d;
    end;

    if ( FLAG__data_phascor2d ),
      meas.data_phascor2d = data_phascor2d;
%      meas.scan_phascor2d = scan_phascor2d;
    end;

    if ( FLAG__offline ),
      meas.offline = offline;
    end;

    if ( FLAG__noiseadjscan ),
      meas.noiseadjscan = noiseadjscan;
    end;

    if ( FLAG__patrefscan ),
      meas.patrefscan = patrefscan;
      clear patrefscan
%      meas.scan_patrefscan = scan_patrefscan;
    end;

    if ( FLAG__patrefscan_phascor ),
      meas.patrefscan_phascor = patrefscan_phascor;
    end;

    if ( FLAG__smsrefscan ),
      meas.smsrefscan = smsrefscan;
    end;

    if ( FLAG__smsrefscan_phascor ),
      meas.smsrefscan_phascor = smsrefscan_phascor;
    end;

    if ( FLAG__phasestabscan ),
      meas.refphasestabscan = refphasestabscan;
      meas.refphasestabtime = refphasestabtime;
      meas.phasestabscan = phasestabscan;
      meas.phasestabtime = phasestabtime;
    end;

    if ( ~isempty(patrefphasestabscan) ),
      meas.patrefphasestabscan = patrefphasestabscan;
      meas.patrefphasestabtime = patrefphasestabtime;
    end;

    if ( DO__STORE_TIMESTAMP ),
      meas.timestamp_data = timestamp_data;
      if ( FLAG__data_phascor1d ),
        meas.timestamp_data_phascor1d = timestamp_data_phascor1d;
      end;
    end;

    meas.lastpos = pos;

    % if DO__MDH_SAVE == 0 this mdh is only the last one, but still can be helpful
    if ( DO__MDH_SAVE > -1 ),
      meas.mdh = mdh;
    end;

    if ( exist('read_meas_prot', 'file') && IS__RELEASE_VERSION && ~DO__EXTRACT_REPETITION ),
      [meas.prot, meas.evp] = read_meas_prot(wildfile, header);

      if ( DO__CANONICAL_REORDER_COIL_CHANNELS ),
        meas.prot.sCoilElementID_tElement = {meas.prot.sCoilElementID_tElement{coil_order}};
      end;


      if ( FLAG__phasestabscan ),
        meas.prot.alTS = TS * 1000;  % integer-valued, so microseconds (same as TE, TR)
      end;
    
      if ( FLAG__allow_asym ),
	meas.evp.MDH_ushSamplesInScan = ushSamplesInScan;
	meas.prot.flAsymmetry_DERIVED = samples_asym;

      end;
    
    end;

    
    if ( DO__HDR_SAVE ),
      write_meas_dat__hdr(measfile, header);
    end;

%    meas.timestamp = timestamp;

    meas.options = options;

    if ( FLAG__status_OK ),
      meas.STATUS = 'success';
    else,
      meas.STATUS = 'ABORTED';
    end;

    meas.data_begin = data_begin;
    meas.head_begin = header_start;

    
    meas.isvd = IS__VD_PLATFORM;
    
    if ( nargout > 0 ),
      varargout{1} = meas;
    end;


  else,

    % if auxiliary data is not present, set arrays to empty arrays in case
    % user requested arrays as output arguments

    if ( ~FLAG__data_phascor1d ),
      data_phascor1d          = [];
    end;

    if ( ~FLAG__data_phascor2d ),
      data_phascor2d          = [];
    end;

    if ( ~FLAG__offline ),
      offline                 = [];
    end;

    if ( ~FLAG__noiseadjscan ),
      noiseadjscan            = [];
    end;

    if ( ~FLAG__patrefscan ),
      patrefscan              = [];
    end;

    if ( ~FLAG__patrefscan_phascor ),
      patrefscan_phascor      = [];
    end;

    if ( ~FLAG__phasestabscan ),
      refphasestabscan       = [];
      refphasestabtime       = [];
      phasestabscan          = [];
      phasestabtime          = [];
    end;


    if ( nargout > 0 ),
      varargout{1} = meas.data;
      varargout{2} = data_phascor1d;
      varargout{3} = data_phascor2d;
      varargout{4} = noiseadjscan;
      varargout{5} = patrefscan;
      varargout{6} = patrefscan_phascor;
      varargout{7} = phasestabscan;
      varargout{8} = refphasestabscan;
    end;

  end;  %% IF DO__RETURN_STRUCT

  return;


%**************************************************************************%
function bit = read_meas__extract_bit(mdh_eval_info, FLAG__EvalInfoMask)

  bit = bitget(mdh_eval_info, FLAG__EvalInfoMask+1);

  return;


%**************************************************************************%
function flag = read_meas__extract_flag(mdh_eval_info, flag_str)

%slow!  persistent EvalInfoMask
  EvalInfoMask = read_meas__evaluation_info_mask_definition;

  field_str = sprintf('MDH_%s', flag_str);
  if ( ~isfield(EvalInfoMask, field_str) ),
    error('flag %s is not a valid EvalInfoMask flag', flag_str);
  end;

  flag = bitget(mdh_eval_info, EvalInfoMask.(field_str)+1);


  return;


%**************************************************************************%
function mdh_struct = read_meas__mdh_struct(varargin)
% snarfed from <n4/pkg/MrServers/MrMeasSrv/SeqIF/MDH/mdh.h>

  mdh_struct = struct(...
      'ulFlagsAndDMALength', [], ...
      'lMeasUID', [], ...
      'ulScanCounter', [], ...
      'ulTimeStamp', [], ...
      'ulPMUTimeStamp', [], ...
      'aulEvalInfoMask', [], ...
      'sEvalInfoMask', [], ...
      'ushSamplesInScan', [], ...
      'ushUsedChannels', [], ...
      'ushLine', [], ...
      'ushAcquisition', [], ...   % note: acquisition is same as average
      'ushSlice', [], ...
      'ushPartition', [], ...
      'ushEcho', [], ...
      'ushPhase', [], ...
      'ushRepetition', [], ...
      'ushSet', [], ...
      'ushSeg', [], ...
      'ushIda', [], ...
      'ushIdb', [], ...
      'ushIdc', [], ...
      'ushIdd', [], ...
      'ushIde', [], ...
      'ushPre', [], ...
      'ushPost', [], ...
      'ushKSpaceCentreColumn', [], ...
      'ushDummy', [], ...
      'fReadOutOffcentre', [], ...
      'ulTimeSinceLastRF', [], ...
      'ushKSpaceCentreLineNo', [], ...
      'ushKSpaceCentrePartitionNo', [], ...
      'aushIceProgramPara', [], ...
      'aushFreePara', [], ...
      'flSag', [], ...
      'flCor', [], ...
      'flTra', [], ...
      'aflQuaternion', [], ...
      'ushChannelId', [], ...
      'ushPTABPosNeg', []);

  return;


%**************************************************************************%
function EvalInfoMask = read_meas__evaluation_info_mask_definition(varargin)
% snarfed from <n4/pkg/MrServers/MrMeasSrv/SeqIF/MDH/MdhProxy.h>
% //##ModelId=3AFAAF7801CF

  EvalInfoMask.MDH_ACQEND            = 0;
  EvalInfoMask.MDH_RTFEEDBACK        = 1;
  EvalInfoMask.MDH_HPFEEDBACK        = 2;
  EvalInfoMask.MDH_ONLINE            = 3;
  EvalInfoMask.MDH_OFFLINE           = 4;
  EvalInfoMask.MDH_SYNCDATA          = 5;   % readout contains synchronous data
  EvalInfoMask.six                   = 6;
  EvalInfoMask.seven                 = 7;
  EvalInfoMask.MDH_LASTSCANINCONCAT  = 8;   % Flag for last scan in concatenation
  EvalInfoMask.nine                  = 9;

  EvalInfoMask.MDH_RAWDATACORRECTION = 10;  % Correct the rawdata with the rawdata correction factor
  EvalInfoMask.MDH_LASTSCANINMEAS    = 11;  % Flag for last scan in measurement
  EvalInfoMask.MDH_SCANSCALEFACTOR   = 12;  % Flag for scan specific additional scale factor
  EvalInfoMask.MDH_2NDHADAMARPULSE   = 13;  % 2nd RF excitation of HADAMAR
  EvalInfoMask.MDH_REFPHASESTABSCAN  = 14;  % reference phase stabilization scan
  EvalInfoMask.MDH_PHASESTABSCAN     = 15;  % phase stabilization scan
  EvalInfoMask.MDH_D3FFT             = 16;  % execute 3D FFT
  EvalInfoMask.MDH_SIGNREV           = 17;  % sign reversal
  EvalInfoMask.MDH_PHASEFFT          = 18;  % execute phase fft
  EvalInfoMask.MDH_SWAPPED           = 19;  % swapped phase/readout direction
  EvalInfoMask.MDH_POSTSHAREDLINE    = 20;  % shared line
  EvalInfoMask.MDH_PHASCOR           = 21;  % phase correction data
  EvalInfoMask.MDH_PATREFSCAN        = 22;  % additional scan for PAT reference line/partition
  EvalInfoMask.MDH_PATREFANDIMASCAN  = 23;  % additional scan for PAT reference line/partition that is also used as image scan
  EvalInfoMask.MDH_REFLECT           = 24;  % reflect line
  EvalInfoMask.MDH_NOISEADJSCAN      = 25;  % noise adjust scan --> Not used in NUM4
  EvalInfoMask.MDH_SHARENOW          = 26;  % all lines are acquired from the actual and previous e.g. phases
  EvalInfoMask.MDH_LASTMEASUREDLINE  = 27;  % indicates that the current line is the last measured line of all succeeding e.g. phases
  EvalInfoMask.MDH_FIRSTSCANINSLICE  = 28;  % indicates first scan in slice (needed for time stamps)
  EvalInfoMask.MDH_LASTSCANINSLICE   = 29;  % indicates last scan in slice (needed for time stamps)
  EvalInfoMask.MDH_TREFFECTIVEBEGIN  = 30;  % indicates the begin time stamp for TReff (triggered measurement)
  EvalInfoMask.MDH_TREFFECTIVEEND    = 31;  % indicates the end time stamp for TReff (triggered measurement)

  EvalInfoMask.thirty_two            = 32;
  EvalInfoMask.thirty_three          = 33;
  EvalInfoMask.thirty_four           = 34;
  EvalInfoMask.thirty_five           = 35;
  EvalInfoMask.thirty_six            = 36;
  EvalInfoMask.thirty_seven          = 37;
  EvalInfoMask.thirty_eight          = 38;
  EvalInfoMask.thirty_nine           = 39;

  EvalInfoMask.MDH_FIRST_SCAN_IN_BLADE       = 40;  % Marks the first line of a blade
  EvalInfoMask.MDH_LAST_SCAN_IN_BLADE        = 41;  % Marks the last line of a blade
  EvalInfoMask.MDH_LAST_BLADE_IN_TR          = 42;  % Set for all lines of the last BLADE in each TR interval


  EvalInfoMask.MDH_RETRO_LASTPHASE           = 45;  % Marks the last phase in a heartbeat
  EvalInfoMask.MDH_RETRO_ENDOFMEAS           = 46;  % Marks an ADC at the end of the measurement
  EvalInfoMask.MDH_RETRO_REPEATTHISHEARTBEAT = 47;  % Repeat the current heartbeat when this bit is found
  EvalInfoMask.MDH_RETRO_REPEATPREVHEARTBEAT = 48;  % Repeat the previous heartbeat when this bit is found
  EvalInfoMask.MDH_RETRO_ABORTSCANNOW        = 49;  % Just abort everything
  EvalInfoMask.MDH_RETRO_LASTHEARTBEAT       = 50;  % This adc is from the last heartbeat (a dummy)
  EvalInfoMask.MDH_RETRO_DUMMYSCAN           = 51;  % This adc is just a dummy scan, throw it away
  EvalInfoMask.MDH_RETRO_ARRDETDISABLED      = 52;  % Disable all arrhythmia detection when this bit is found


  return;


  % VB13A ICE manual

  %% There are 16 ICE Dimensions listed in the table below:
  %% Dimension   Description
  %% COL         Column related to pixel property ICE_PP_FIX
  %% LIN         Line related to pixel property ICE_PP_FIX
  %% CHA         Element related to pixel property ICE_PP_FIX
  %% SET         Set related to pixel property ICE_PP_FIX
  %% ECO         Contrast (echo) related to pixel property ICE_PP_VAR
  %% PHS         Phase related to pixel property ICE_PP_FIX
  %% REP         Measurement Repeat (Repetition) related to pixel property ICE_PP_VAR
  %% SEG         Segment related to pixel property ICE_PP_FIX
  %% PAR         Partition related to pixel property ICE_PP_VAR
  %% SLC         Slice related to pixel property ICE_PP_VAR
  %% IDA         1st free dispose of an ICE Dimension related to pixel property ICE_PP_VAR
  %% IDB         2nd free dispose of an ICE Dimension related to pixel property ICE_PP_VAR
  %% IDC         3rd free dispose of an ICE Dimension related to pixel property ICE_PP_VAR
  %% IDD         4th free dispose of an ICE Dimension related to pixel property ICE_PP_FIX
  %% IDE         5th free dispose of an ICE Dimension related to pixel property ICE_PP_FIX
  %% AVE         Average related to pixel property ICE_PP_VAR
  %%
  %% Table 14: 16 ICE Dimensions in the order of generation in the IceProF.
  %% There are 11 IceDimension with an intended use (e.g. SLC for slices)
  %% and further 5 IceDimensions which can be used freely.


%**************************************************************************%
function read_meas__dump_flags(mdh_eval_info)
% for debugging:  read_meas__dump_flags(aulEvalInfoMask)

  EvalInfoMask = read_meas__evaluation_info_mask_definition;

  flag_array = fieldnames(EvalInfoMask);

  % stop at 32 since EvalInfoMask is cast as a uint32
  for ind = 1:32,

    flag_str = flag_array{ind};

    if ( isempty(strfind(flag_str, 'MDH_')) ),
      continue;
    end;

    flag_str = regexprep(flag_str, 'MDH_', '');

    if ( read_meas__extract_flag(mdh_eval_info, flag_str) )
      disp(flag_str);
    end;

  end;

  return;


%**************************************************************************%
function disp_optional(str)
% locally-defined version of "disp" that can be silenced with a global variable


  global QUIET;   if ( isempty(QUIET) ), QUIET = 0; end;

  if ( QUIET ),
    return;
  end;

  disp(str);


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/read_meas_dat.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
