function  Merkmale = harris_detektor(Image, varargin) 
    % In dieser Funktion soll der Harris-Detektor implementiert werden, der
    % Merkmalspunkte aus dem Bild extrahiert
    
    % check whether the picture does only contain one colour dimension
    if size(Image,3) ~= 1
        return;
    end

    % initial values 
    segment_length=3;
    k=0.05;
    tau=1;
    do_plot = false;

    % read optional arguments
    n=1;
    while n+1<=nargin
        switch varargin{n}
            case 'do_plot'
                do_plot = varargin{n+1};
            case 'segment_length'
                segment_length = varargin{n+1};
                if mod(segemnt_length,2) == 0
                    disp 'Error: segment length is an even number! we should never got here';
                    return;
                end
            case 'k'
                k = varargin{n+1};
            case 'tau'
                tau = varargin{n+1};
            otherwise
                disp 'Error: unknown input in harris detector! we should never got here'
        end
        n= n+2;
    end
    fprintf('segment_length = %d, k = %d, tau = %d, do_plot = %d\n',segment_length,k,tau,do_plot);

    % allocate memory
    Merkmale = zeros(size(Image));
    H = zeros(size(Image));

    % compute gradient
    [Fx,Fy] = sobel_xy(Image);

    % run through image
    disp 'compute G ...'
    p=(segment_length-1)/2;
    y=1+p;
    while(y<size(Image,2)-p)
        x=1+p;
        while(x<size(Image,1)-p)
            % Compute the approximated Harris-Matrix for the Pixel Image(x,y)
            G = zeros(2);
            for ix = -p:p
                for iy = -p:p
                    w = weight(ix, iy, segment_length);
                    G = G + w*[Fx(x+ix,y+iy),Fy(x+ix,y+iy)]'*[Fx(x+ix,y+iy),Fy(x+ix,y+iy)];
                end
            end
            H(x,y) = det(G)-k*(trace(G)^2);
            
            if H(x,y)>tau
                Merkmale(x,y)=1;% Ecke (black)
            %elseif H(x,y)<-tau
                %Merkmale(x,y)=2;% Kante
            %else
                %Merkmale(x,y)=3;% Flaeche
            end
            x=x+1;
        end
        y=y+1;
    end
    disp 'computation of G finished!';
    
    % plot features
    if (do_plot==true)
        plot_harris(Merkmale, Image);
        
    end
end

function w = weight(ix,iy,segment_length)
    % Computes the Weight, such that pixels that are closer to the the
    % central pixel (i.e. ix and iy are small) have a higher weight
    if (segment_length == 3)
    switch abs(ix*iy)
        case 0
            w = 4;
        case 1
            w = 2;
        case 2
            w = 1;
        otherwise
            disp('error');
    end
    end
end

function plot_harris(Merkmale, Image)
%% -------------------------------------
    % Plots the image
    figure;
    imshow(Merkmale);

end