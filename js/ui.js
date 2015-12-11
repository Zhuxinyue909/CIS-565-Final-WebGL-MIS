var cfg;

(function() {
    'use strict';

    var Cfg = function() {
        // TODO: Define config fields and defaults here
        this.debugView = -1;
        this.intensity=50;
        this.GammaCorrection=true;
        this.pos_x=-0.5;
        this.pos_y=0;
        this.pos_z=-3
        this.enableEffect0 = false;
    };

    var init = function() {
        cfg = new Cfg();

        var gui = new dat.GUI();
        // TODO: Define any other possible config values
        gui.add(cfg, 'debugView', {
            'None':             -1,
            '0 Depth':           0,
            '1 Position':        1,
            '2 Geometry normal': 2,
            '3 Color map':       3,
            '4 Normal map':      4,
            '5 Surface normal':  5
        });
       // gui.add(cfg, 'debugScissor');

        var eff0 = gui.addFolder('EFFECT NAME HERE');
        eff0.open();
        eff0.add(cfg, 'intensity').min(0).max(100).step(1);
        eff0.add(cfg, 'GammaCorrection');
        eff0.add(cfg,'pos_x').min(-1.0).max(4).step(0.2);
        eff0.add(cfg,'pos_y').min(0).max(4).step(0.2);
        eff0.add(cfg,'pos_z').min(-4).max(4).step(0.2);
        // TODO: add more effects toggles and parameters here
    };

    window.handle_load.push(init);
})();
