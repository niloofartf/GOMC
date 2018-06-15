
class t_state
{
    public:
        //! Constructor
        t_state();
  
        int                        natoms;         //!< Number of atoms, local + non-local; this is the size of \p x, \p v and \p cg_p, when used
        std::vector<std::vector<float>> x;              //!< The coordinates (natoms)
        // Center of Mass
        // Box Dimensions
        // Cell List

};

