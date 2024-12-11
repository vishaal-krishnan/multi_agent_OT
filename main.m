% n_samples = 500;
% N = 100;
% n_neigh = 10;
% n_trial_samples = 300;
% 
% c = 0.05;
% T = 400;
% 
% %---------------Target distribution------------------
% dom_size = 5;
% x = (-dom_size:dom_size/100:dom_size);
% y = (-dom_size:dom_size/100:dom_size);
% [X,Y] = meshgrid(x,y);
% components = 20;
% mean_des = 8*(rand(components,2) - 0.5);
% var_des = [0.25,0.25];
% gm = gmdistribution(mean_des, var_des);
% samples = random(gm,n_samples);
% 
% sizegrid = size(X);
% Z = zeros(size(X));
% 
% for i=1:sizegrid(1)
%     for j=1:sizegrid(2)
%         Z(i,j) = pdf(gm, [X(i,j), Y(i,j)]) ;
%     end
% end
% 
% for i=1:sizegrid(1)
%     for j=1:sizegrid(2)
%         Z(i,j) = pdf(gm, [X(i,j), Y(i,j)]) ;
%     end
% end
% 
% mean_init = [5, 5];
% var_init = [2,6];
% gm_init = gmdistribution(mean_init, var_init);
% pos_init = random(gm_init,N);

n_steps = 1;
count = zeros(N,T,n_steps);
% t = (1:1:n_steps);
% t = [1,10,20,30,40,50];
t = 20;

for s=1:n_steps
    
    D = 0;
    phi = zeros(N,1);
    pos = pos_init;
    pos_save = zeros(T,N,2);
    
    for i=1:T
        
        pos_save(i,:,:) = pos;
        
        id_samples_pos = knnsearch(pos, samples);
        count(:,i,s) = (histcounts(id_samples_pos,(0.1:1:N+0.1)')/n_samples)';
        
        for j=1:t
            id_pos_pos = knnsearch(pos,pos,'K',n_neigh);
            A = zeros(N,N);
            for k=1:N
                A(k,id_pos_pos(k,2:n_neigh)) = 1;
            end
            G = graph(0.5*(A+A'));
            L = laplacian(G);
            
            phi = phi - (1/(n_neigh+1))*L*phi + (1/N)*ones(N,1) - count(:,i);
        end
        
        PHI = scatteredInterpolant(pos(:,1),pos(:,2),phi(:,1),'linear');
        
        mean_trial = [0, 0];
        var_trial = [0.05,0.05];
        gm_trial = gmdistribution(mean_trial, var_trial);
        samples_trial = random(gm_trial,n_trial_samples);
        
        for m=1:N
            d = (samples_trial(:,1).*samples_trial(:,1) + samples_trial(:,2).*samples_trial(:,2)).^(0.5);
            
            samples_geq_trial = samples_trial((d>=0.1),:);
            %         samples_geq_trial = samples_trial(:,:);
            samples_geq_size = size(samples_geq_trial);
            
            pos_trial = zeros(samples_geq_size(1),2);
            
            pos_trial(:,1) = samples_geq_trial(:,1) + pos(m,1);
            pos_trial(:,2) = samples_geq_trial(:,2) + pos(m,2);
            
            
            [temp, ind] = min(d((d>=0.1)) + PHI(pos_trial(:,1),pos_trial(:,2)));
            %         [temp, ind] = min(PHI(pos_trial(:,1),pos_trial(:,2)));
            pos_temp = pos_trial(ind,:);
            %        pos_temp = min(pos(m,:),pos_trial(ind,:));
            D = D + ((pos_temp - pos(m,:))*(pos_temp - pos(m,:))')^(0.5);
            %        pos(m,:) = min(pos(m,:),pos_trial(ind,:));
            pos(m,:) = pos_trial(ind,:);
            
        end
        
    end
    
end



%--------------------------------PLOTS-------------------------------------

hold on;

for s=1:n_steps
    
    % pts = linspace(-2*dom_size, 2*dom_size, 50);
    % N = histcounts2(samples(:,2), samples(:,1), pts, pts);
    % colormap(flipud(gray))
    % imagesc(pts, pts, N);
    % hold on;
    % scatter(pos_save(100,:,1), pos_save(100,:,2),'r+','LineWidth',2.5);
    % axis equal;
    % set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]));
    % grid minor;
    
    
    plot((1:1:T),var(count(:,:,s)),'LineWidth',1.5);
    xlabel('T', 'FontSize', 18);
    ylabel('Var($\mu^*(\mathcal{V}_i)$)','Interpreter','latex','FontSize', 18);
    xlim([0,300]);
    ylim([0, max(var(count(:,:,s)))]);
    grid on;
    
end

var_count_20 = [var_count_20,var(count(:,:,s))'];

% for i=1:T
%     
%     clf;
%     
%     subplot(1,2,1)
%     
% %     surf(X,Y,-Z);
% %     colormap(gray);
% %     shading interp;
% %     view(2);
% %     hold on;
% %      
% %     [posx, posy] = voronoi(pos_save(i,:,1),pos_save(i,:,2));
% %     plot(pos_save(i,:,1), pos_save(i,:,2), 'r+','LineWidth',1.5,'MarkerSize',5);    
% %     axis(2*dom_size*[-1 1 -1 1]);
% %     grid off;
% 
%     pts = linspace(-2*dom_size, 2*dom_size, 50);
%     N = histcounts2(samples(:,2), samples(:,1), pts, pts);
%     colormap(flipud(gray))
%     imagesc(pts, pts, N);
%     hold on; 
%     scatter(pos_save(i,:,1), pos_save(i,:,2),'r+','LineWidth',2.5);
%     axis equal;
%     set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]));
%     grid minor;
%  
%     hold off;
%     
%     subplot(1,2,2)
%     plot((0:1:i-1)/10,var_time_count(1:i),'LineWidth',1.5);
%     xlabel('T', 'FontSize', 18);
%     ylabel('Var($\mu^*(\mathcal{V}_i)$)','Interpreter','latex','FontSize', 18);
%     xlim([0,40]);
%     ylim([0, max(var_time_count)]);
%     
%     F(i) = getframe(gcf);
%     drawnow
%       
% end
%         
% %%%% 1fps
% writerObj = VideoWriter('testvideo.avi');
% writerObj.FrameRate = 10;
% 
% %open video writer
% open(writerObj);
% 
% for i=1:length(F)
%     frame = F(i);
%     writeVideo(writerObj, frame);
% end
% close(writerObj);
