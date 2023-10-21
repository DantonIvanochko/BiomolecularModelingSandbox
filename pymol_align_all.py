from pymol import cmd

tmp_object_list = cmd.get_names()
for i in range(1,len(tmp_object_list)):
    cmd.align(tmp_object_list[i], tmp_object_list[0])
