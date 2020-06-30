function muscles=whichMuscles(muscle)

switch muscle
    case 1 % triceps
        muscles=[129:133];%[104:108 129:133];
    case 2 % deltoids
        muscles=[37:51];
    case 3 % lat dorsi
        muscles=[90:95];
    case 4 % serratus anterior
        muscles=[25:36];
    case 5 % biceps/brachialis
        muscles=[83:85];
    case 6 % supra/infraspinatus
        muscles=[55:60 68:71];
    case 7 % rhomboids
        muscles=[20:24];
    case 8 % lower pectoralis
        muscles=[96:101];
    case 9 % upper pectoralis
        muscles=[102:103];
    case 10 % static
        muscles=[];
end

