function  Merkmale = harris_detektor(Image, varargin) 
    % In dieser Funktion soll der Harris-Detektor implementiert werden, der
    % Merkmalspunkte aus dem Bild extrahiert
    
    % check whether the picture does only contain one colour dimension
    if size(Image,3) ~= 1
        return;
    end

    % initial values 
    segment_length=2;
    k=0.05;
    tau=1;
    do_plot = false;

    % read optional arguments
    n=1;
    while n<=nargin
        switch varargin{n}
            case 'do_plot'
                do_plot = varargin{n+1};
            case 'segment_length'
                segment_length = varargin{n+1};
            case 'k'
                k = varargin{n+1};
            case 'tau'
                tau = varargin{n+1};
            otherwise
                disp 'Error: unknown input in harris detector! we should never got here'
        end
        n= n+2;
    end

    disp(segment_length, k, tau, do_plot);

    % allocate memory
    Merkmale = zeros(size(Image));

    % compute gradient
    [Fx,Fy] = sobel_xy(Image);

    % run through image
    disp 'compute G ...'
    y=-(segment_length-1)/2;
    while(y<size(Image,2)-(segment_length-1)/2)
        x=-(segment_length-1)/2;
        while(x<size(Image,1)-(segment_length-1)/2)
            % Berechne die approximierte Harris Matrix G (steht fuer
            % Aenderung des Bildsegments)
            G = zeros(2,2);
            for ix = -(segment_length-1)/2:(segment_length-1)/2
                for iy = -(segment_length-1)/2:(segment_length-1)/2
                    w = weight(ix, iy, segment_length);
                    G = G + w*[Fx(x+ix,y+iy),Fy(x+ix,y+iy)]'*[Fx(x+ix,y+iy),Fy(x+ix,y+iy)];
                end
            end
            H = det(G)-k*(tr(G))^2;
            if H>tau
                Merkmale(x,y)='side';
            elseif H<tau
                Merkmale(x,y)='angle';
            else
                Merkmale(x,y)='area';
            end
            x=x+segment_length;
        end
        y=y+segment_length;
    end
    disp 'computing of G finished!'
    
    % plot features
    if (do_plot==true)
        plot_harris(Merkmale, Image);
    end
end

function w = weight(ix,iy,segment_length)
    % Computes the Weight, such that pixels that are closer to the the
    % central pixel (i.e. ix and iy are small) have a higher weight
    w = 1;
end

function plot_harris(Merkmale, Image)
    % Plots the image
end