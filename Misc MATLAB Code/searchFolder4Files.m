function files = searchFolder4Files(origin)

folder = dir(origin);
files  = '';

   for n = 3:size(folder,1)
       suborigin = [origin,filesep,folder(n).name];
       switch isdir(suborigin)
        case 1
            files = getintofolder(suborigin,files);
        case 0
            files{size(files,1)+1,1} = suborigin;
       end
   end
end

   function files = getintofolder(origin,files)

       folder = dir(origin);
       for n = 3:size(folder,1)
           suborigin = [origin,filesep,folder(n).name];
           switch isdir(suborigin)
            case 1
                files = getintofolder(suborigin,files);
            case 0
                files{size(files,1)+1,1} = suborigin;
           end
       end
   end