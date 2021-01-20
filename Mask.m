classdef Mask
    %MASK Object that describes a masking operation on a (set of) frames
    %   This object combines code for packing/unpacking frames into column
    %   vectors (see pack, unpack), code for generating coordinate ranges
    %   (properties x and y), and cropping frames (crop).
    %
    %   Currently, only binary masks are supported.
    %   TODO: allow apodized (smooth) masks.
    %
    properties
        x_range
        y_range
        bounds
        x %coordinates from -1 to 1
        y %coordinates from -1 to 1
        bitmask %binary mask
        filter %smoothing filter
        count %number of non-zero elements in the mask
    end
    
    methods
        function obj = Mask(bounds, varargin)
            %MASK Construct a data mask that can be used to crop, pack and unpack data
            %bounds =   four element vector of form [y x h w] (top-left
            %           position, width and height in pixel coordinates)
            %           for the rectangular region where the mask is placed.
            %
            % remaining options are given as name-value pairs
            % 'Shape'   = Mask type
            %           'circle' (default), gives a circular (or elliptical) mask
            %           'ring', a ring (see InnerRadius)
            % 'InnerRadius' = Only used when shape is 'ring'. Specifies
            %               the excluded part) of the ring, this value
            %               should be between 0 (full disk is used) and 1
            %               (ring has zero width: no pixels selected). 
            %
            % 'Filter'  = Filter to smooth edges of the mask. 
            %           'cosine' (default) = uses cosine (tukey) window
            %
            % 'Smoothness' = Parameter specifying how much the edges are
            %           smoothed. 0 (default) = hard edges.
            %
            
            %% Parse and validate attributes
            p = inputParser;
            p.addRequired('bounds', @(v)validateattributes(v, {'numeric'}, {'vector', 'numel', 4}));            
            p.addParameter('Shape', 'circular', @ischar);
            p.addParameter('InnerRadius', 0.5, @isnumeric); %used for ring shape only
            p.addParameter('Filter', 'cosine', @ischar);
            p.addParameter('Smoothness', 0, @isnumeric);
            p.parse(bounds, varargin{:});
            opt = p.Results;
            
            %% set up coordinates
            obj.x = linspace(-1, 1, bounds(4));
            obj.y = linspace(-1, 1, bounds(3)).';
            obj.x_range = bounds(2) + (1:bounds(4));
            obj.y_range = bounds(1) + (1:bounds(3));
            obj.bounds = bounds;

            %% construct mask
            r2 = (obj.x.^2 + obj.y.^2); %radius squared
            if strcmp(opt.Shape, 'circular')
                obj.bitmask = r2 <= 1;
            elseif strcmp(opt.Shape, 'ring')
                obj.bitmask = (r2 <= 1) & (r2 > opt.InnerRadius^2);
            else
                error('Unsupported mask type. Must be ''circular'' or ''ring''.');
            end
            obj.count = sum(obj.bitmask(:) ~= 0);
            
            %% Optionally apply a smoothing filter to the outside of the mask
            if opt.Smoothness > 0
                s = min(opt.Smoothness, 1);
                r = min(max((sqrt(r2) + (s-1)) / s, 0), 1);
                if strcmp(opt.Filter, 'cosine')
                    f = 0.5 + 0.5 * cos(r * pi);
                else
                    error('Unsupported window type. Must be ''cosine''.');
                end
            else
                f = 1;
            end
            obj.filter = obj.bitmask .* f;
        end
        
        function Y = pack(obj, varargin)
            % Convert array data to packed format. 
            % PACK(x) 
            %
            % x is a single frame or a multi-dimensional array holding
            % multiple frames. The first two dimensions of x correspond to
            % the size of a single frame.
            %
            % When the frame size is larger than the mask size, the
            % frame is cropped with the bounds given in the Mask constructor.
            % When the frame size equals the mask size, it is assumed that
            % the data was already cropped.
            %
            % The function discards all values of the frame for which mask=0 and 
            % then stores each frame as a column in the 2-D output array y.
            %
            % PACK(x1, x2, ...)
            % Packs multiple arrays together. This is useful, for example
            % for combining x and y gradients into a single matrix.
            % equivalent to cat(2, mask.pack(x1), mask.pack(x2), ...)
            %
            X = varargin{1};
            sz = size(X);
            N_frames = prod(sz(3:end));
            
            if nargin > 2 %handle multi-input case
                Ys = cell(nargin-1, 1);
                for j=1:(nargin-1)
                    X = varargin{j};
                    if any(size(X) ~= sz)
                        error('When packing multiple arrays, sizes must match.');
                    end
                    Ys{j} = pack(obj, X);
                end
                Y = vertcat(Ys{:});
                return
            end
            
            % crop data if needed
            X = crop(obj, X);
            
            % pack each frame into a single column of y
            Y = zeros(obj.count, N_frames, 'like', X);
            X = reshape(X, [numel(X)/N_frames, N_frames]);
            for n=1:size(Y,2)
                Y(:, n) = X(obj.bitmask(:), n);
            end
        end
        
        function varargout = unpack(obj, Y)
            % unpacks previously packed data. The result is a 3-D array
            % where each slice corresponds to a frame.
            % If multiple arrays were packed together (see pack), each
            % of them is returned as a separate argument
            %
            if mod(size(Y, 1), obj.count)  ~= 0
                error('Incompatible size');
            end
            
            % Special case for unpacking multiple arrays
            if size(Y, 1) > obj.count
                varargout{nargout} = []; %pre-allocate output cell array
                for j=1:nargout % max(nargout, size(y, 1) / obj.count)
                    varargout{j} = unpack(obj, Y((j-1) * obj.count + (1:obj.count), :));
                end
                return
            end
            
            % Unpack single frames one at a time
            N_frames = size(Y,2);
            X = zeros(numel(obj.bitmask), N_frames, 'like', Y);
            for n=1:N_frames
                X(obj.bitmask(:), n) = Y(:,n);
            end
            varargout{1} = reshape(X, [size(obj.bitmask), N_frames]);
        end
        
        function X = crop(obj, X)
            % Crops matrix to the size of the bounding box of the mask.
            % If data is already cropped, this function does nothing
            % If the data is smaller than the bounding box, an error 
            % occurs.
            %
            if size(X, 2) ~= obj.bounds(4) || size(X, 1) ~= obj.bounds(3)
                X = X(obj.y_range, obj.x_range);
            end
        end
        
        function X = apply(obj, X)
            % Applies the mask to frame or set of frames X
            X = crop(obj, X) .* obj.filter;
        end
    end
end

