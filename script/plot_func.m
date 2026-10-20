

function Y = plot_func(dimension,matrix,x, y, z,  fontsize, linewidth, x_dis, y_dis, width, height, x_label, y_label, title1)

%% 绘制单个图形
if matrix(1) == 0
    if dimension == 2
        figure('color',[1 1 1], 'position', [x_dis, y_dis, width, height]); % [位置横纵 宽 高]
        plot([0 0],[0 0],'r','LineWidth',1.5);hold on;
        set(gca,'FontSize',fontsize,'FontName','Times New Roman', 'LineWidth', linewidth);
        xlabel(x_label,'FontName','Times New Roman','FontSize',fontsize);
        ylabel(y_label,'FontName','Times New Roman','FontSize',fontsize);
        title(title1,'FontName','Monospaced','FontSize',fontsize);
        axis tight;
        for ix = 1:length(y(1,:))
            plot(x, y(:,ix),  'Color', rand(1,3), 'LineWidth', linewidth);hold on;
        end
        xlim([min(x), max(x)]);
        
    end
    if dimension == 3
        figure('color',[1 1 1], 'position', [x_dis, y_dis, width, height]); % [位置横纵 宽 高]
        h = pcolor(x, y, z); % 绘制图像，并返回句柄
        set(h, 'EdgeColor', 'none');
        axis tight; % 自动调整坐标轴的范围，使其紧贴着数据的最小值和最大值
        colorbar;
        colormap jet;
        set(gca, 'TickDir', 'in', 'Box', 'off', 'LineWidth', linewidth);
        set(gca,'FontName','Times New Roman','FontSize',fontsize); % Times New Roman
        xlabel(x_label,'FontName','Times New Roman','FontSize',fontsize);
        ylabel(y_label,'FontName','Times New Roman','FontSize',fontsize);
        title(title1,'FontName','Monospaced','FontSize',fontsize);
    end
end


%% 一图绘制多个图形
if matrix(1) ~= 0
    if matrix(4) == 1
        figure('color',[1 1 1], 'position', [x_dis, y_dis, width, height]);
    end
    if dimension == 2
    subplot(matrix(1),matrix(2),matrix(3));
        plot([0 0],[0 0],'r','LineWidth',1.5);hold on;
        set(gca,'FontSize',fontsize,'FontName','Times New Roman', 'LineWidth', linewidth);
        xlabel(x_label,'FontName','Times New Roman','FontSize',fontsize);
        ylabel(y_label,'FontName','Times New Roman','FontSize',fontsize);
        title(title1,'FontName','Monospaced','FontSize',fontsize);
        axis tight;
        for ix = 1:length(y(1,:))
            plot(x, y(:,ix),  'Color', rand(1,3), 'LineWidth', linewidth);hold on;
        end
        xlim([min(x), max(x)]);hold on;
    end

    if dimension == 3
        subplot(matrix(1),matrix(2),matrix(3));
%     figure('color',[1 1 1], 'position', [x_dis, y_dis, width, height]); % [位置横纵 宽 高]
    h = pcolor(x, y, z); % 绘制图像，并返回句柄
    set(h, 'EdgeColor', 'none');
    axis tight; % 自动调整坐标轴的范围，使其紧贴着数据的最小值和最大值
%     colorbar;
    colormap jet;
    set(gca, 'TickDir', 'in', 'Box', 'off', 'LineWidth', linewidth);
    set(gca,'FontName','Times New Roman','FontSize',fontsize); % Times New Roman
    xlabel(x_label,'FontName','Times New Roman','FontSize',fontsize);
    ylabel(y_label,'FontName','Times New Roman','FontSize',fontsize);
    title(title1,'FontName','Monospaced','FontSize',fontsize);
    end
end


end










