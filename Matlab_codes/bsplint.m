function [idx, val] = bsplint(X, Y)

%     [IDX, VAL] = BSPLINT(X, Y) Returns the knots 'IDX' and values 'VAL'
%     of the intersections between the Bézier curves defined from the set
%     of control points 'X' and 'Y'. The following cases are supported:
% 
%     - BSPLINT(X, Y)    Intersections between curves 'X' and 'Y'
%     - BSPLINT(X)       Self-intersections of 'X'
%     - BSPLINT(X, X)    Self-intersections of 'X'
%     - BSPLINT(X, n)    Intersections between 'X' and a scalar 'N'


    % Parse input arguments
    [a, c] = size(X);
    if nargin == 1
        if a < 4
            idx = [];
            val = [];
            return
        end
        Y = X;
    end
    
    % Adjust intersecting curve
    [b, d] = size(Y);
    
    % Special case: scalar intersection
    if b == 1 && d == 1
        t = (0 : 1 / (a - 1) : 1)';
        idx = zeros((a - 1) * c, 1);
        val = 0;
        
        % Intersect each dimension
        for i = 1 : c
            x = X(:, i);
            ii = bsplint([t, x], [0, Y; 1, Y]);
            k = size(ii, 1);
            if k > 0
                idx(val + 1 : val + k) = ii(:, 1);
                val = val + k;
            end
        end
        
        % Evaluate intersections
        idx = unique(idx(idx > 0));
        if isempty(idx)
            val = [];
            return
        end
        if nargout > 1
            k = binomial(a);
            val = bspl(X, idx, a, k);
        end
        return
    end
    if c ~= d
        if c > d
            Y = [Y, zeros(b, c - d)];
        else
            X = [X, zeros(a, d - c)];
            c = d;
        end
    end
    
    % One-Dimensional curve case
    if c == 1
        X = [(0 : 1 / (a - 1) : 1)', X];
        Y = [(0 : 1 / (b - 1) : 1)', Y];
        c = 2;
    end
    si = isequal(X, Y);
    
    % Initialise parameters
    t = 1e-12;
    n = 10 * (a - 1) * (b - 1) / (1 + si);
    n = min(max(n, 200), 1e5);
    
    % Initialise knots
    u = 0 : 1 / (n - 1) : 1;
    v = u;
    
    % Initialise binomial coefficients
    f = binomial(a);
    g = binomial(b);
    
    % Get overlapping bounding boxes
    x = bspl(X, u, a, f);
    y = bspl(Y, v, b, g);
    k = bbint(x, y, n, c);
    
    % De-cluster adjacent cases
    k = decls(k);
    if si
        
        % Self-intersections: remove mirrors
        k = issym(k);
    end
    
    % Return if no overapping regions
    idx = [];
    m = size(k, 1);
    if m == 0
        val = [];
        return
    end
    
    % Loop over the number of occurrences
    for ii = m : -1 : 1
        d = k(ii, 1);
        e = k(ii, 2);
        
        % Avoid same knot span in self-intersections
        if si
            if abs(d - e) < 5
                continue
            end
        end
        
        % Filter goodness of intersection
        r = x(d, :) - y(e, :);
        r = r .* r;
        r = sum(r, 2);
        if r < t
            idx(ii, :) = [u(d), v(e)];
            
            % Special case: intersection at origin
            if d * e == 1
                idx(ii, :) = t;
            end
            continue
        end
        
        % Refine query knot span
        w = 100;
        k(ii, :) = [u(d), v(e)];
        p = kspan(d, 5, u, w);
        q = kspan(e, 5, v, w);
        
        % Loop until fit tolerance is met
        i = 0;
        z = false;
        while r > t
            
            % Overlap bounding boxes
            xx = bspl(X, p, a, f);
            yy = bspl(Y, q, b, g);
            j = bbint(xx, yy, w, c);
            if isempty(j)
                break
            end
            
            % Deculster data
            j = decls(j);
            if si
                j = issym(j);
            end
            
            % Correct duplicates
            if size(j, 1) > 1
                if z
                    idx(ii, :) = 0;
                    break
                else
                    z = true;
                    j = fix(sum(j) / size(j, 1));
                end
            end
            
            % Refine knot span if residuals decrease
            s = xx(j(1), :) - yy(j(2), :);
            s = s .* s;
            s = sum(s, 2);
            
            % Update enhanced values
            if s < r
                r = s;
                k(ii, :) = [p(j(1)), q(j(2))];
                idx(ii, :) = k(ii, :);
                p = kspan(j(1), 2, p, w);
                q = kspan(j(2), 2, q, w);
                
            % Refine knot resolution otherwise
            else
                s = w;
                w = fix(1.2 * w);
                p = p(1) : (p(s) - p(1)) / (w - 1) : p(s);
                q = q(1) : (q(s) - q(1)) / (w - 1) : q(s);
            end
            
            % Exit condition
            i = i + 1;
            if i == 5
                break
            end
        end
    end
    
    % Parse indices
    idx = idx(any(idx, 2), :);
    if isempty(idx)
        val = [];
        return
    end
    
    % Remove duplicates
    if size(idx, 1) > 1
        if si

            % Check linear indices
            for i = size(idx, 1) : -1 : 1
                j = idx(i, :) - idx;
                j = j .* j;
                j = sum(j, 2);
                j(i, :) = 1;
                j = find(j < 1e-4);
                if any(j)
                    idx(j, :) = NaN;

                % Check mirrowed values
                else
                    j = idx(i, [2, 1]) - idx;
                    j = j .* j;
                    j = sum(j, 2);
                    j(i, :) = 1;
                    j = find(j < 1e-4);
                    if any(j)
                        idx(i, :) = NaN;
                    end
                end
            end
            idx(any(isnan(idx), 2), :) = [];

        % Remove close indices
        else
            i = diff(idx);
            i = i .* i;
            i = sum(i, 2);
            idx(i < 1e-4, :) = [];
        end
    else
        if isempty(idx)
            val = [];
            return
        end
    end
    
    % Parse outputs
    if nargout > 1
        val = bspl(X, idx(:, 1), a, f);
    end
end

% Obtain binomial coefficients
function b = binomial(n)

    % Table look-up for linear coefficients
    if n < 59
        b = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,3,3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,4,6,4,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,5,10,10,5,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,6,15,20,15,6,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,7,21,35,35,21,7,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,8,28,56,70,56,28,8,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,9,36,84,126,126,84,36,9,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,10,45,120,210,252,210,120,45,10,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,11,55,165,330,462,462,330,165,55,11,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,12,66,220,495,792,924,792,495,220,66,12,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,13,78,286,715,1287,1716,1716,1287,715,286,78,13,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,14,91,364,1001,2002,3003,3432,3003,2002,1001,364,91,14,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,15,105,455,1365,3003,5005,6435,6435,5005,3003,1365,455,105,15,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,16,120,560,1820,4368,8008,11440,12870,11440,8008,4368,1820,560,120,16,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,17,136,680,2380,6188,12376,19448,24310,24310,19448,12376,6188,2380,680,136,17,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,18,153,816,3060,8568,18564,31824,43758,48620,43758,31824,18564,8568,3060,816,153,18,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,19,171,969,3876,11628,27132,50388,75582,92378,92378,75582,50388,27132,11628,3876,969,171,19,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,20,190,1140,4845,15504,38760,77520,125970,167960,184756,167960,125970,77520,38760,15504,4845,1140,190,20,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,21,210,1330,5985,20349,54264,116280,203490,293930,352716,352716,293930,203490,116280,54264,20349,5985,1330,210,21,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,22,231,1540,7315,26334,74613,170544,319770,497420,646646,705432,646646,497420,319770,170544,74613,26334,7315,1540,231,22,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,23,253,1771,8855,33649,100947,245157,490314,817190,1144066,1352078,1352078,1144066,817190,490314,245157,100947,33649,8855,1771,253,23,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,24,276,2024,10626,42504,134596,346104,735471,1307504,1961256,2496144,2704156,2496144,1961256,1307504,735471,346104,134596,42504,10626,2024,276,24,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,25,300,2300,12650,53130,177100,480700,1081575,2042975,3268760,4457400,5200300,5200300,4457400,3268760,2042975,1081575,480700,177100,53130,12650,2300,300,25,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,26,325,2600,14950,65780,230230,657800,1562275,3124550,5311735,7726160,9657700,10400600,9657700,7726160,5311735,3124550,1562275,657800,230230,65780,14950,2600,325,26,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,27,351,2925,17550,80730,296010,888030,2220075,4686825,8436285,13037895,17383860,20058300,20058300,17383860,13037895,8436285,4686825,2220075,888030,296010,80730,17550,2925,351,27,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,28,378,3276,20475,98280,376740,1184040,3108105,6906900,13123110,21474180,30421755,37442160,40116600,37442160,30421755,21474180,13123110,6906900,3108105,1184040,376740,98280,20475,3276,378,28,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,29,406,3654,23751,118755,475020,1560780,4292145,10015005,20030010,34597290,51895935,67863915,77558760,77558760,67863915,51895935,34597290,20030010,10015005,4292145,1560780,475020,118755,23751,3654,406,29,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,30,435,4060,27405,142506,593775,2035800,5852925,14307150,30045015,54627300,86493225,119759850,145422675,155117520,145422675,119759850,86493225,54627300,30045015,14307150,5852925,2035800,593775,142506,27405,4060,435,30,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,31,465,4495,31465,169911,736281,2629575,7888725,20160075,44352165,84672315,141120525,206253075,265182525,300540195,300540195,265182525,206253075,141120525,84672315,44352165,20160075,7888725,2629575,736281,169911,31465,4495,465,31,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,32,496,4960,35960,201376,906192,3365856,10518300,28048800,64512240,129024480,225792840,347373600,471435600,565722720,601080390,565722720,471435600,347373600,225792840,129024480,64512240,28048800,10518300,3365856,906192,201376,35960,4960,496,32,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,33,528,5456,40920,237336,1107568,4272048,13884156,38567100,92561040,193536720,354817320,573166440,818809200,1037158320.00000,1166803110.00000,1166803110.00000,1037158320.00000,818809200,573166440,354817320,193536720,92561040,38567100,13884156,4272048,1107568,237336,40920,5456,528,33,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,34,561,5984,46376,278256,1344904,5379616,18156204,52451256,131128140,286097760,548354040,927983760,1391975640.00000,1855967520.00000,2203961430.00000,2333606220.00000,2203961430.00000,1855967520.00000,1391975640.00000,927983760,548354040,286097760,131128140,52451256,18156204,5379616,1344904,278256,46376,5984,561,34,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,35,595,6545,52360,324632,1623160,6724520,23535820,70607460,183579396,417225900,834451800,1476337800.00000,2319959400.00000,3247943160.00000,4059928950.00000,4537567650.00000,4537567650.00000,4059928950.00000,3247943160.00000,2319959400.00000,1476337800.00000,834451800,417225900,183579396,70607460,23535820,6724520,1623160,324632,52360,6545,595,35,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,36,630,7140,58905,376992,1947792,8347680,30260340,94143280,254186856,600805296,1251677700.00000,2310789600.00000,3796297200.00000,5567902560.00000,7307872110.00000,8597496600.00000,9075135300.00000,8597496600.00000,7307872110.00000,5567902560.00000,3796297200.00000,2310789600.00000,1251677700.00000,600805296,254186856,94143280,30260340,8347680,1947792,376992,58905,7140,630,36,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,37,666,7770,66045,435897,2324784,10295472,38608020,124403620,348330136,854992152,1852482996.00000,3562467300.00000,6107086800.00000,9364199760.00000,12875774670.0000,15905368710.0000,17672631900.0000,17672631900.0000,15905368710.0000,12875774670.0000,9364199760.00000,6107086800.00000,3562467300.00000,1852482996.00000,854992152,348330136,124403620,38608020,10295472,2324784,435897,66045,7770,666,37,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,38,703,8436,73815,501942,2760681,12620256,48903492,163011640,472733756,1203322288.00000,2707475148.00000,5414950296.00000,9669554100.00000,15471286560.0000,22239974430.0000,28781143380.0000,33578000610.0000,35345263800.0000,33578000610.0000,28781143380.0000,22239974430.0000,15471286560.0000,9669554100.00000,5414950296.00000,2707475148.00000,1203322288.00000,472733756,163011640,48903492,12620256,2760681,501942,73815,8436,703,38,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,39,741,9139,82251,575757,3262623,15380937,61523748,211915132,635745396,1676056044.00000,3910797436.00000,8122425444.00000,15084504396.0000,25140840660.0000,37711260990.0000,51021117810.0000,62359143990.0000,68923264410.0000,68923264410.0000,62359143990.0000,51021117810.0000,37711260990.0000,25140840660.0000,15084504396.0000,8122425444.00000,3910797436.00000,1676056044.00000,635745396,211915132,61523748,15380937,3262623,575757,82251,9139,741,39,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,40,780,9880,91390,658008,3838380,18643560,76904685,273438880,847660528,2311801440.00000,5586853480.00000,12033222880.0000,23206929840.0000,40225345056.0000,62852101650.0000,88732378800.0000,113380261800.000,131282408400.000,137846528820.000,131282408400.000,113380261800.000,88732378800.0000,62852101650.0000,40225345056.0000,23206929840.0000,12033222880.0000,5586853480.00000,2311801440.00000,847660528,273438880,76904685,18643560,3838380,658008,91390,9880,780,40,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,41,820,10660,101270,749398,4496388,22481940,95548245,350343565,1121099408.00000,3159461968.00000,7898654920.00000,17620076360.0000,35240152720.0000,63432274896.0000,103077446706.000,151584480450.000,202112640600.000,244662670200.000,269128937220.000,269128937220.000,244662670200.000,202112640600.000,151584480450.000,103077446706.000,63432274896.0000,35240152720.0000,17620076360.0000,7898654920.00000,3159461968.00000,1121099408.00000,350343565,95548245,22481940,4496388,749398,101270,10660,820,41,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,42,861,11480,111930,850668,5245786,26978328,118030185,445891810,1471442973.00000,4280561376.00000,11058116888.0000,25518731280.0000,52860229080.0000,98672427616.0000,166509721602.000,254661927156.000,353697121050.000,446775310800.000,513791607420.000,538257874440.000,513791607420.000,446775310800.000,353697121050.000,254661927156.000,166509721602.000,98672427616.0000,52860229080.0000,25518731280.0000,11058116888.0000,4280561376.00000,1471442973.00000,445891810,118030185,26978328,5245786,850668,111930,11480,861,42,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,43,903,12341,123410,962598,6096454,32224114,145008513,563921995,1917334783.00000,5752004349.00000,15338678264.0000,36576848168.0000,78378960360.0000,151532656696.000,265182149218.000,421171648758.000,608359048206.000,800472431850.000,960566918220.000,1052049481860.00,1052049481860.00,960566918220.000,800472431850.000,608359048206.000,421171648758.000,265182149218.000,151532656696.000,78378960360.0000,36576848168.0000,15338678264.0000,5752004349.00000,1917334783.00000,563921995,145008513,32224114,6096454,962598,123410,12341,903,43,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,44,946,13244,135751,1086008,7059052,38320568,177232627,708930508,2481256778.00000,7669339132.00000,21090682613.0000,51915526432.0000,114955808528.000,229911617056.000,416714805914.000,686353797976.000,1029530696964.00,1408831480056.00,1761039350070.00,2012616400080.00,2104098963720.00,2012616400080.00,1761039350070.00,1408831480056.00,1029530696964.00,686353797976.000,416714805914.000,229911617056.000,114955808528.000,51915526432.0000,21090682613.0000,7669339132.00000,2481256778.00000,708930508,177232627,38320568,7059052,1086008,135751,13244,946,44,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,45,990,14190,148995,1221759,8145060,45379620,215553195,886163135,3190187286.00000,10150595910.0000,28760021745.0000,73006209045.0000,166871334960.000,344867425584.000,646626422970.000,1103068603890.00,1715884494940.00,2438362177020.00,3169870830126.00,3773655750150.00,4116715363800.00,4116715363800.00,3773655750150.00,3169870830126.00,2438362177020.00,1715884494940.00,1103068603890.00,646626422970.000,344867425584.000,166871334960.000,73006209045.0000,28760021745.0000,10150595910.0000,3190187286.00000,886163135,215553195,45379620,8145060,1221759,148995,14190,990,45,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,46,1035,15180,163185,1370754,9366819,53524680,260932815,1101716330.00000,4076350421.00000,13340783196.0000,38910617655.0000,101766230790.000,239877544005.000,511738760544.000,991493848554.000,1749695026860.00,2818953098830.00,4154246671960.00,5608233007146.00,6943526580276.00,7890371113950.00,8233430727600.00,7890371113950.00,6943526580276.00,5608233007146.00,4154246671960.00,2818953098830.00,1749695026860.00,991493848554.000,511738760544.000,239877544005.000,101766230790.000,38910617655.0000,13340783196.0000,4076350421.00000,1101716330.00000,260932815,53524680,9366819,1370754,163185,15180,1035,46,1,0,0,0,0,0,0,0,0,0,0,0,0,0;1,47,1081,16215,178365,1533939,10737573,62891499,314457495,1362649145.00000,5178066751.00000,17417133617.0000,52251400851.0000,140676848445.000,341643774795.000,751616304549.000,1503232609098.00,2741188875414.00,4568648125690.00,6973199770790.00,9762479679106.00,12551759587422.0,14833897694226.0,16123801841550.0,16123801841550.0,14833897694226.0,12551759587422.0,9762479679106.00,6973199770790.00,4568648125690.00,2741188875414.00,1503232609098.00,751616304549.000,341643774795.000,140676848445.000,52251400851.0000,17417133617.0000,5178066751.00000,1362649145.00000,314457495,62891499,10737573,1533939,178365,16215,1081,47,1,0,0,0,0,0,0,0,0,0,0,0,0;1,48,1128,17296,194580,1712304,12271512,73629072,377348994,1677106640.00000,6540715896.00000,22595200368.0000,69668534468.0000,192928249296.000,482320623240.000,1093260079344.00,2254848913647.00,4244421484512.00,7309837001104.00,11541847896480.0,16735679449896.0,22314239266528.0,27385657281648.0,30957699535776.0,32247603683100.0,30957699535776.0,27385657281648.0,22314239266528.0,16735679449896.0,11541847896480.0,7309837001104.00,4244421484512.00,2254848913647.00,1093260079344.00,482320623240.000,192928249296.000,69668534468.0000,22595200368.0000,6540715896.00000,1677106640.00000,377348994,73629072,12271512,1712304,194580,17296,1128,48,1,0,0,0,0,0,0,0,0,0,0,0;1,49,1176,18424,211876,1906884,13983816,85900584,450978066,2054455634.00000,8217822536.00000,29135916264.0000,92263734836.0000,262596783764.000,675248872536.000,1575580702584.00,3348108992991.00,6499270398159.00,11554258485616.0,18851684897584.0,28277527346376.0,39049918716424.0,49699896548176.0,58343356817424.0,63205303218876.0,63205303218876.0,58343356817424.0,49699896548176.0,39049918716424.0,28277527346376.0,18851684897584.0,11554258485616.0,6499270398159.00,3348108992991.00,1575580702584.00,675248872536.000,262596783764.000,92263734836.0000,29135916264.0000,8217822536.00000,2054455634.00000,450978066,85900584,13983816,1906884,211876,18424,1176,49,1,0,0,0,0,0,0,0,0,0,0;1,50,1225,19600,230300,2118760,15890700,99884400,536878650,2505433700.00000,10272278170.0000,37353738800.0000,121399651100.000,354860518600.000,937845656300.000,2250829575120.00,4923689695575.00,9847379391150.00,18053528883775.0,30405943383200.0,47129212243960.0,67327446062800.0,88749815264600.0,108043253365600,121548660036300,126410606437752,121548660036300,108043253365600,88749815264600.0,67327446062800.0,47129212243960.0,30405943383200.0,18053528883775.0,9847379391150.00,4923689695575.00,2250829575120.00,937845656300.000,354860518600.000,121399651100.000,37353738800.0000,10272278170.0000,2505433700.00000,536878650,99884400,15890700,2118760,230300,19600,1225,50,1,0,0,0,0,0,0,0,0,0;1,51,1275,20825,249900,2349060,18009460,115775100,636763050,3042312350.00000,12777711870.0000,47626016970.0000,158753389900.000,476260169700.000,1292706174900.00,3188675231420.00,7174519270695.00,14771069086725.0,27900908274925.0,48459472266975.0,77535155627160.0,114456658306760,156077261327400,196793068630200,229591913401900,247959266474052,247959266474052,229591913401900,196793068630200,156077261327400,114456658306760,77535155627160.0,48459472266975.0,27900908274925.0,14771069086725.0,7174519270695.00,3188675231420.00,1292706174900.00,476260169700.000,158753389900.000,47626016970.0000,12777711870.0000,3042312350.00000,636763050,115775100,18009460,2349060,249900,20825,1275,51,1,0,0,0,0,0,0,0,0;1,52,1326,22100,270725,2598960,20358520,133784560,752538150,3679075400.00000,15820024220.0000,60403728840.0000,206379406870.000,635013559600.000,1768966344600.00,4481381406320.00,10363194502115.0,21945588357420.0,42671977361650.0,76360380541900.0,125994627894135,191991813933920,270533919634160,352870329957600,426384982032100,477551179875952,495918532948104,477551179875952,426384982032100,352870329957600,270533919634160,191991813933920,125994627894135,76360380541900.0,42671977361650.0,21945588357420.0,10363194502115.0,4481381406320.00,1768966344600.00,635013559600.000,206379406870.000,60403728840.0000,15820024220.0000,3679075400.00000,752538150,133784560,20358520,2598960,270725,22100,1326,52,1,0,0,0,0,0,0,0;1,53,1378,23426,292825,2869685,22957480,154143080,886322710,4431613550.00000,19499099620.0000,76223753060.0000,266783135710.000,841392966470.000,2403979904200.00,6250347750920.00,14844575908435.0,32308782859535.0,64617565719070.0,119032357903550,202355008436035,317986441828055,462525733568080,623404249591760,779255311989700,903936161908052,973469712824056,973469712824056,903936161908052,779255311989700,623404249591760,462525733568080,317986441828055,202355008436035,119032357903550,64617565719070.0,32308782859535.0,14844575908435.0,6250347750920.00,2403979904200.00,841392966470.000,266783135710.000,76223753060.0000,19499099620.0000,4431613550.00000,886322710,154143080,22957480,2869685,292825,23426,1378,53,1,0,0,0,0,0,0;1,54,1431,24804,316251,3162510,25827165,177100560,1040465790.00000,5317936260.00000,23930713170.0000,95722852680.0000,343006888770.000,1108176102180.00,3245372870670.00,8654327655120.00,21094923659355.0,47153358767970.0,96926348578605.0,183649923622620,321387366339585,520341450264090,780512175396135,1.08592998315984e+15,1.40265956158146e+15,1.68319147389775e+15,1.87740587473211e+15,1.94693942564811e+15,1.87740587473211e+15,1.68319147389775e+15,1.40265956158146e+15,1.08592998315984e+15,780512175396135,520341450264090,321387366339585,183649923622620,96926348578605.0,47153358767970.0,21094923659355.0,8654327655120.00,3245372870670.00,1108176102180.00,343006888770.000,95722852680.0000,23930713170.0000,5317936260.00000,1040465790.00000,177100560,25827165,3162510,316251,24804,1431,54,1,0,0,0,0,0;1,55,1485,26235,341055,3478761,28989675,202927725,1217566350.00000,6358402050.00000,29248649430.0000,119653565850.000,438729741450.000,1451182990950.00,4353548972850.00,11899700525790.0,29749251314475.0,68248282427325.0,144079707346575,280576272201225,505037289962205,841728816603675,1.30085362566023e+15,1.86644215855598e+15,2.48858954474130e+15,3.08585103547921e+15,3.56059734862986e+15,3.82434530038022e+15,3.82434530038022e+15,3.56059734862986e+15,3.08585103547921e+15,2.48858954474130e+15,1.86644215855598e+15,1.30085362566023e+15,841728816603675,505037289962205,280576272201225,144079707346575,68248282427325.0,29749251314475.0,11899700525790.0,4353548972850.00,1451182990950.00,438729741450.000,119653565850.000,29248649430.0000,6358402050.00000,1217566350.00000,202927725,28989675,3478761,341055,26235,1485,55,1,0,0,0,0;1,56,1540,27720,367290,3819816,32468436,231917400,1420494075.00000,7575968400.00000,35607051480.0000,148902215280.000,558383307300.000,1889912732400.00,5804731963800.00,16253249498640.0,41648951840265.0,97997533741800.0,212327989773900,424655979547800,785613562163430,1.34676610656588e+15,2.14258244226390e+15,3.16729578421620e+15,4.35503170329728e+15,5.57444058022051e+15,6.64644838410907e+15,7.38494264901008e+15,7.64869060076044e+15,7.38494264901008e+15,6.64644838410907e+15,5.57444058022051e+15,4.35503170329728e+15,3.16729578421620e+15,2.14258244226390e+15,1.34676610656588e+15,785613562163430,424655979547800,212327989773900,97997533741800.0,41648951840265.0,16253249498640.0,5804731963800.00,1889912732400.00,558383307300.000,148902215280.000,35607051480.0000,7575968400.00000,1420494075.00000,231917400,32468436,3819816,367290,27720,1540,56,1,0,0,0;1,57,1596,29260,395010,4187106,36288252,264385836,1652411475.00000,8996462475.00000,43183019880.0000,184509266760.000,707285522580.000,2448296039700.00,7694644696200.00,22057981462440.0,57902201338905.0,139646485582065,310325523515700,636983969321700,1.21026954171123e+15,2.13237966872931e+15,3.48934854882978e+15,5.30987822648010e+15,7.52232748751348e+15,9.92947228351779e+15,1.22208889643296e+16,1.40313910331192e+16,1.50336332497705e+16,1.50336332497705e+16,1.40313910331192e+16,1.22208889643296e+16,9.92947228351779e+15,7.52232748751348e+15,5.30987822648010e+15,3.48934854882978e+15,2.13237966872931e+15,1.21026954171123e+15,636983969321700,310325523515700,139646485582065,57902201338905.0,22057981462440.0,7694644696200.00,2448296039700.00,707285522580.000,184509266760.000,43183019880.0000,8996462475.00000,1652411475.00000,264385836,36288252,4187106,395010,29260,1596,57,1,0,0;1,58,1653,30856,424270,4582116,40475358,300674088,1916797311.00000,10648873950.0000,52179482355.0000,227692286640.000,891794789340.000,3155581562280.00,10142940735900.0,29752626158640.0,79960182801345.0,197548686920970,449972009097765,947309492837400,1.84725351103293e+15,3.34264921044054e+15,5.62172821755909e+15,8.79922677530988e+15,1.28322057139936e+16,1.74517997710313e+16,2.21503612478474e+16,2.62522799974487e+16,2.90650242828897e+16,3.00672664995410e+16,2.90650242828897e+16,2.62522799974487e+16,2.21503612478474e+16,1.74517997710313e+16,1.28322057139936e+16,8.79922677530988e+15,5.62172821755909e+15,3.34264921044054e+15,1.84725351103293e+15,947309492837400,449972009097765,197548686920970,79960182801345.0,29752626158640.0,10142940735900.0,3155581562280.00,891794789340.000,227692286640.000,52179482355.0000,10648873950.0000,1916797311.00000,300674088,40475358,4582116,424270,30856,1653,58,1,0;1,59,1711,32509,455126,5006386,45057474,341149446,2217471399.00000,12565671261.0000,62828356305.0000,279871768995.000,1119487075980.00,4047376351620.00,13298522298180.0,39895566894540.0,109712808959985,277508869722315,647520696018735,1.39728150193517e+15,2.79456300387033e+15,5.18990272147347e+15,8.96437742799963e+15,1.44209549928690e+16,2.16314324893035e+16,3.02840054850248e+16,3.96021610188786e+16,4.84026412452961e+16,5.53173042803384e+16,5.91322907824307e+16,5.91322907824307e+16,5.53173042803384e+16,4.84026412452961e+16,3.96021610188786e+16,3.02840054850248e+16,2.16314324893035e+16,1.44209549928690e+16,8.96437742799963e+15,5.18990272147347e+15,2.79456300387033e+15,1.39728150193517e+15,647520696018735,277508869722315,109712808959985,39895566894540.0,13298522298180.0,4047376351620.00,1119487075980.00,279871768995.000,62828356305.0000,12565671261.0000,2217471399.00000,341149446,45057474,5006386,455126,32509,1711,59,1];
        b = b(n, 1 : n);
        
    % Central terms for logarithmic approach
    else
        b = gammaln(2 : n - 1) + gammaln(n - 1 : -1 : 2);
        b = gammaln(n) - b';
    end
end

% Generate Bézier curve
function y = bspl(x, t, n, b)
   
    % Low order: Bernstein matrix
    if n < 59
        t = repmat(t(:), 1, n);
        n = n - 1;
        y = t .^ (0 : n);
        y = y .* (1 - t) .^ (n : -1 : 0);
        y = y .* b;
        
        % Evaluate curve
        y = y * x;
        
    % High order: Logarithmic construction
    else
        
        % Initialise variables
        t = t(:);
        n = n - 1;
        m = length(t);
        a = abs(min(x(:))) + 1;
        x = x + a;
        
        % Initialise Bézier curve
        y = x(1, :) .* (1 - t) .^ n;
        y = y + x(n + 1, :) .*  t .^ n;
        
        % Initialise logarithmic terms
        x(1, :) = [];
        x(n, :) = [];
        x = log(x) + b;
        
        % Obtain logarithmic Bernstein matrix
        n = n - 1;
        t = log(t) .* (1 : n) + log(1 - t) .* (n : -1 : 1);
        
        % Evaluate control points
        t = repmat(t', size(x, 2), 1);
        t = t + x(:);
        
        % Convert to linear form
        t = exp(t);
        t = reshape(t', m, n, []);
        t = permute(t, [1, 3, 2]);
        t = reshape(t, [], n);
        t = sum(t, 2);
        t = reshape(t, m, []);
        
        % Evaluate curve
        y = y + t;
        y = y - a;
    end
end

% Overlap bounding boxes
function k = bbint(x, y, k, n)

    % Get linear indicies
    i = 1 : k - 1;
    j = 2 : k;
    
    % Opposite vertex superposition in each dimension
    k = min(x(i, 1), x(j, 1)) <= max(y(i, 1), y(j, 1))';
    k = k & max(x(i, 1), x(j, 1)) >= min(y(i, 1), y(j, 1))';
    for n = 2 : n
        k = k & min(x(i, n), x(j, n)) <= max(y(i, n), y(j, n))';
        k = k & max(x(i, n), x(j, n)) >= min(y(i, n), y(j, n))';
    end
    
    % Get overlapping indices
    [i, j] = find(k);
    k = [i, j];
end

% Decluster intersecting indicies
function k = decls(k)
    if size(k, 1) < 2
        return
    end
    
    % Mark to delete adjacent boxes
    i = sortrows(k);
    k = abs(diff(i));
    k = sum(k, 2);
    k = fix(0.5 * k) > 2;
    k = find([true; k]);
    k = [i(k, 1), i(k, 2)];
end

% Remove mirrors
function k = issym(k)
    if size(k, 1) < 2
        return
    end
    
    % Remove symmetric knots and decluster
    i = abs(k(:, 1) - k(:, 2)) > 2;
    k = decls(k(i, :));
    
    % Loop for symmetries in flipped indices
    for s = size(k, 1) : -1 : 1
        i = k(:, [2, 1]);
        i = abs(k(s, :) - i);
        i = sum(i, 2);
        i = fix(i / s) < 2;
        if any(i)
            k(s, :) = 0;
        end
    end
    
    % Keep unique values
    k = k(any(k, 2), :);
end

% Refine knot span
function u = kspan(k, n, u, m)
    k = k + n * [-1, 1];
    n = length(u);
    k(k < 1) = 1;
    k(k > n) = n;
    u = u(k(1)) : (u(k(2)) - u(k(1))) / (m - 1) : u(k(2));
end