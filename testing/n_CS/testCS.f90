Program testCS

Use Kinds, Only: dp
Use FileIO_utilities, Only: slash
Use n_Cross_Sections, Only:  CS_type
Use n_Cross_Sections, Only:  Setup_Cross_Sections

Implicit None

Type(CS_type) :: CS

CS = Setup_Cross_Sections( & 
                           & resources_directory = '', & 
                           & cs_setup_file = '', & 
                           & elastic_only = .FALSE., & 
                           & aniso_dist = .TRUE., & 
                           & E_min = 0._dp , E_max = 30000._dp)

End Program testCS
