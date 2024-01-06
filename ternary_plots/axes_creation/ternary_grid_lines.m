function handle = ternary_grid_lines( handle, varargin )
%ternary_grid_lines plots a set of grid lines for all three axes
%   
%   grid_pnts is structure array grid_pnts(1:3).values(1:N), where N
%   can vary for each axis. Ticklinelength gives any extra length on the
%   base axis used for ticks
%   
    
    %% Check Inputs
    
    % Check if handle is given
    if (nargin==0)
       error('Too few inputs') 
    end
    
    % Create defaul settings if no extra arguments supplied
    if ( nargin < 2 || isempty( varargin ) )
        varargin = {'LineStyle','-','LineWidth',1.5,'Color',[0 0 0 0.4]};
    end
    
    %% Loop Axes and Set grid lines
    
    % Extract a copy of grid points from handle
    try
        grid_pnts = handle.grid.grid_pnts;
    catch
       error('Spacing array, handle.grid.grid_pnts, not defined'); 
    end
    
    % Local copy of wlimits 
    wlimits = handle.grid.wlimits;
    
    % Loop each axis
    for iaxis = 1:3
        
        % Loop grid lines
        for i = 1:numel( grid_pnts(iaxis).values )    
            
            % First point at the base, second on the far side, both stored in A/B
            [E,F,~] = tern2base( iaxis, grid_pnts(iaxis).values(i), ...
                                       wlimits, 0.0 );
            
            % Plot3 based on two points in A/B, pass varargin elements
            handle.grid.lines(i,iaxis) = ternary_plot3( wlimits, 1, E, 2, F, [], varargin{:} );
            
        end
        
        % Link properties related to formatting to each axes
        props = {'LineStyle','LineWidth','Visible','Marker',...
           'MarkerEdgeColor','MarkerFaceColor','MarkerIndices',...
           'Selected','MarkerSize'};
            
        % Link these proprties to all Gridlines on an axis
        handle.grid.link_lines(iaxis) = linkprop(handle.grid.lines(:,iaxis),props);
        
    end
    
    % Link All gridlines to Z position for each of changing later
     handle.grid.link_axes = linkprop(handle.grid.lines(:),'ZData');
    
end
