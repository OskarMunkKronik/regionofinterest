clc
clear variables;close all force
% Create a parfor wait bar

% Number of Steps 
NbrePts = 50;

% Waitbar's Message Field
Msg = 'Progress...!';

% Create ParFor Waitbar
[hWaitbar,hWaitbarMsgQueue]= ParForWaitbarCreateMH(Msg,NbrePts);

parfor i1=1:NbrePts
    % Send data from worker to client
    hWaitbarMsgQueue.send(0);
    pause(0.1)
end
delete(hWaitbarMsgQueue);
close(hWaitbar)