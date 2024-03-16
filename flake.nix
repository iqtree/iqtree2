{
  description = "IQ-TREE 2";

  inputs.flake-utils.url = "github:numtide/flake-utils";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable-small";

  outputs =
    { self
    , flake-utils
    , nixpkgs
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        # # For Nixpkgs and the NUR repository, the following two forms are
        # # equivalent because the flake just imports the base directory.
        pkgs = nixpkgs.legacyPackages.${system};
        # pkgs = import nixpkgs {
        #   inherit system;
        #   config = { allowUnfree = true; };
        # };
      in
      {
        devShell = pkgs.mkShell {
          packages = with pkgs; [
            # See https://github.com/NixOS/nixpkgs/issues/59209.
            bashInteractive
          ];
          nativeBuildInputs = with pkgs; [ cmake ];
          buildInputs = with pkgs; [ boost eigen zlib ];

          cmakeFlags = [
            "-DIQTREE_FLAGS=omp"
          ];
        };
      }
    );
}
