function pconfig = cbm_config(d,algorithm,pconfig)

if nargin<3 || isempty(pconfig)
    pconfig = struct('d',d);
end

p = inputParser;
p.addParamValue('algorithm',algorithm);
p.addParamValue('d',d);
p.addParamValue('verbose',1);
p.addParamValue('functionname','',@(arg)ischar(arg)|| isempty(arg) );
p.addParamValue('fname_prog',sprintf('cbm_%s_%0.4f.mat',algorithm,now),@valid_fname);
p.addParamValue('flog',1,@valid_flog);
p.addParamValue('save_data',0,@(arg)(arg==1)||(arg==0));
% p.addParamValue('loop',0,@(arg)(arg==1)||(arg==0)); % not used anymore, replaced with algorithm (mandatory)

if any(strcmp({'lap'},algorithm))  
    p.addParamValue('inits',[],@(arg)(ismatrix(arg) && (size(arg,2)==d) ));    
    p.addParamValue('numinit',min(7*d,100),@(arg)isscalar(arg));
    p.addParamValue('save_prog',0,@(arg)(arg==1)||(arg==0));
    p.addParamValue('discard_bad',0,@(arg)isscalar(arg));    
    p.addParamValue('numinit_med',100,@(arg)isscalar(arg));
    p.addParamValue('numinit_up',1000,@(arg)isscalar(arg));        
end
if any(strcmp({'lap','hierlap','loophierlap'},algorithm))
    p.addParamValue('rng',[-5*ones(1,d);5*ones(1,d)],@(arg)valid_rng(d,arg));
    p.addParamValue('largescale','off');
    p.addParamValue('gradient','off');
    p.addParamValue('hessian','off'); 
    p.addParamValue('tolgrad',.001001,@(arg)isscalar(arg));
    p.addParamValue('tolgrad_liberal',0.10,@(arg)(isvector(arg) || isempty(arg)));
    p.addParamValue('reject_Fc',[],@(arg)valid_reject_Fc(arg));
    p.addParamValue('free_group',true(d,1),@(arg)numel(arg)==d);
    p.addParamValue('free_groupvar',true(d,1),@(arg)numel(arg)==d);
end

if any(strcmp({'hierlap','loophierlap'},algorithm))
    p.addParamValue('inits',[],@(arg)(ismatrix(arg) && size(arg,2)==d) || isempty(arg));
    p.addParamValue('numinit',0,@(arg)isscalar(arg));   
    p.addParamValue('numinit_med',10,@(arg)isscalar(arg));
    p.addParamValue('numinit_up',50,@(arg)isscalar(arg));    
    p.addParamValue('tolx',0.05,@(arg)isscalar(arg)); %increase to have faster fitting
    p.addParamValue('tolL',-log(.5),@(arg)isscalar(arg));    
    p.addParamValue('maxiter',50,@(arg)isscalar(arg)); % maximum number of iterations  
    p.addParamValue('terminate','dx',@(arg)valid_terminate(arg));  % termination criteria (default based on mean parameters)
    p.addParamValue('save_prog',1);

end
if strcmp('loophierlap',algorithm)
    p.addParamValue('loop_runtime',30,@(arg)(isscalar(arg)&& (arg<300)));
    p.addParamValue('loop_maxruntime',90,@(arg)(isscalar(arg)&& (arg<600)));
    p.addParamValue('loop_pausesec',30,@(arg)(isscalar(arg)&& (arg<60)));
    p.addParamValue('loop_maxnumrun',3,@(arg)(floor(arg)==arg));    
    p.addParamValue('loop_discard_bad',true,@(arg)(islogical(arg)));
end

if any(strcmp({'slicesample'},algorithm))
    p.addParamValue('inits',[],@(arg)(ismatrix(arg) && size(arg,2)==d));    
    p.addParamValue('nsamp',10^4);
    p.addParamValue('save_prog',1);
    p.addParamValue('burnin',0); % matlab default
    p.addParamValue('thin',1); % matlab default
    p.addParamValue('width',10); % matlab default
end

p.parse(pconfig);
pconfig    = p.Results;

if pconfig.save_prog==0, pconfig.fname_prog = []; end
end

%-------------------------------------------------------
function valid = valid_fname(arg)
valid =1;
if isempty(arg)
    return;
end
try
    save(arg,'valid');
    delete(arg);
catch msg
    warning(msg.message);
    valid = 0;
end

end

function valid = valid_reject_Fc(arg)
valid = 1;
if isempty(arg), return; end
try
    arg(10);
catch %#ok<CTCH>
    valid = 0;
end
end

function valid = valid_terminate(arg)
valid = any(strcmp(arg,{'dL','dx','dxdL'}));
end

function valid = valid_rng(d,arg)
valid = size(arg,1)==2 && (size(arg,2)==d || size(arg,2)==1);
end

function valid = valid_flog(arg)
valid = 0;
if ischar(arg)
    try
        fid = fopen(arg,'w');
        fprintf(fid,'test');
        fclose(fid);
        delete(arg);
        valid = 1;
    catch msg
        warning(msg.message);
        valid = 0;        
    end
    
elseif arg==1
    valid = 1;    
end
end