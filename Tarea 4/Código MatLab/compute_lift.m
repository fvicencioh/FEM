function CL = compute_lift(NODES, TRI, PRESS, INTERIORB, chord)

    npanels = size(INTERIORB, 1);
    CL = 0;

    for k = 1:npanels

        p1 = INTERIORB(k,1);
        p2 = INTERIORB(k,2);

        x1 = NODES(p1,1);  y1 = NODES(p1,2);
        x2 = NODES(p2,1);  y2 = NODES(p2,2);

        % Panel length
        s = norm([x2-x1, y2-y1]);

        % Tangent
        t = ([x2-x1, y2-y1]) / s;

        % Normal (rotación antihoraria)
        n = [ -t(2), t(1) ];

        % Encontrar elemento que contiene ambos nodos
        elem = find(sum(ismember(TRI, [p1 p2]),2) == 2);
        if isempty(elem)
            continue;
        end

        Cp = PRESS(elem);

        % Acumulación de CL
        CL = CL - Cp * n(2) * s / chord;
    end
end