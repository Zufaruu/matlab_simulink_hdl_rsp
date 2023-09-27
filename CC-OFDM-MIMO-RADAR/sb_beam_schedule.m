% Function for making the matrix of sub band allocation per beam and per symbol
% Input variables are the number of the available sub band (Nsb) and the
% number of sub band per beam (Nsb_beam)
% Transmission scheme 1 use Nsb = 8 and Nsb_beam = 3
% Transmission scheme 2 use Nsb = 8 and Nsb_beam = 5
% The condition is that for there are always at least one gap of sub bands
% in every successive beams and every successive symbols
function [sb63_sym]= sb_beam_schedule(Nc,Ncb,Nb,Nsb,Nsb_beam);
sb0 = 1:Nsb;    % Sub band index
sb1 = reshape(reshape(sb0,2,[]).',1,[]);            % Scramble the sub band index for 8 successive beams 
sb63 = reshape(repmat(sb1.',(Nb+1)/Nsb,1),1,[]).';  % Reuse the sb1 combination for the next 8 successive beams
for isb_beam = 1:Nsb_beam
    if isb_beam<5
        sb63_sym (:,isb_beam)= circshift(sb63,-2*(isb_beam-1));
    else
        sb63_sym (:,isb_beam)= circshift(sb63,-1*(isb_beam-4));
    end
end    
