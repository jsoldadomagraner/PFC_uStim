% Dominant and memory subspaces cartoons

savefig = false;

[~, ~, figurepath] = addpaths;

% target location
x1 = [0,1,0,1];
y1 = [0,0,1,1];
z1 = [1,1,1,1]*0.5;

x2 = [0.2,0.8,0.2,0.8];
y2 = [0.2,0.2,0.8,0.8];
z2 = [1,1,1,1]*1.1;

x = x1; y = y1; z = z1;

u = x2-x1; v= y2-y1; w = z2-z1;

[colang, ~,~,~] = plottingspecs;
angorder = [3 4 2 1]; % [225 315 135 45]
colang   = colang(angorder,:);

cols     = lines(5);
colplane = {cols(5,:),cols(4,:)};

s      = 0.12;
alpha  = 1;
alphap = 0.1;

figure;hold on

for ang=1:4
    
    % 3D targets
    [X,Y,Z] = ellipsoid(x(ang),y(ang),z(ang),s,s,s);
    surf(X,Y,Z,'FaceColor',colang(ang,:), ...
        'EdgeColor','none');
    
    xlim([-0.4 1.4])
    ylim([-0.4 1.4])
    zlim([-2.2 1.2])

    XL = get(gca, 'XLim');
    YL = get(gca, 'YLim');
    ZL = get(gca, 'ZLim');
        
    % dominant subspace
    zshift = [0 0 0 0];
    if ang==1
        dom = patch([XL(1), XL(2), XL(2), XL(1)], [YL(1), YL(1), YL(2), YL(2)], ...
        zshift, colplane{1},'FaceAlpha',alphap, ...
        'EdgeColor', colplane{1},'LineWidth',2);

        rotvec    = [0 1 0];    % rotate data around [x,y,z] axis
        rotang    = 50;         % angle of plane in degrees
        rotorigin = [0 0 -0.4]; % point of rotation

        rotate(dom,rotvec,rotang,rotorigin);
    end

    % targets projection
    [X,Y,Z] = ellipsoid(x(ang),y(ang),0,s,s,0);
    domelipses = surf(X,Y,Z,'FaceColor',colang(ang,:),'FAceAlpha',alpha, ...
        'EdgeColor','none');
    
    rotate(domelipses,rotvec,rotang,rotorigin);
   
    % memory subspace
    if ang==1
        zloc = -2;
        fill3([XL(1), XL(2), XL(2), XL(1)], [YL(1), YL(1), YL(2), YL(2)], ...
            [zloc zloc zloc zloc], colplane{2},'FaceAlpha',alphap,...
            'EdgeColor', colplane{2},'LineWidth',2);
    end
    
    % targets projection
    [X,Y,Z] = ellipsoid(x(ang),y(ang),zloc,s,s,0);
    surf(X,Y,Z,'FaceColor',colang(ang,:),'FAceAlpha',alpha, ...
        'EdgeColor','none');
end

ax = gca;
ax.XTick = []; ax.YTick = []; ax.ZTick = [];
rotation = [33.6 17.4];
view(ax,rotation);

pb = pbaspect;
da = daspect;
daspect([1 1 1])


if savefig
    saveas(gcf,[figurepath 'subspaces/' 'subspaces3Dcartoon'],'svg')
end

