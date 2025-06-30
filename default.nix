{ pkgs ? import <nixpkgs> {} }: 
(pkgs.buildFHSEnv {
  name = "sysid";
  targetPkgs = pkgs: (with pkgs; [
    julia
  ]);
  runScript = "bash";
  profile = ''
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${pkgs.julia}/lib/julia
  '';
}).env