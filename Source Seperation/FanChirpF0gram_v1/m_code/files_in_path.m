%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function archivos = files_in_path(ruta,extension,visualizacion)
% FILES_IN_PATH	Obtiene el nombre de los archivos en un directorio.
%	FILES_IN_PATH(RUTA, EXTENSION, VISUALIZACION) Retorna un arreglo de estructuras conteniendo 
%	en la propiedad .name el nombre completo (ruta + nombre de archivo) de los archivos en la ruta
%	indicada por el String RUTA con extensión indicada por el String EXTENSION. La búsqueda es recursiva.
listapath(1).name = ruta;
indicepath = 1;
contpath = dir(ruta);

if(isunix)
    barra = '/';
else
    barra = '\';
end
if(nargin<3)
    visualizacion = 0;
end
    
while indicepath <= length(listapath)
   contpath = dir(listapath(indicepath).name);
   for k = 1:length(contpath)
      if((contpath(k).isdir == 1) & ~strcmp(contpath(k).name,'.') & ~strcmp(contpath(k).name,'..') )
         contpath(k).name;
         listapath(end + 1).name = [listapath(indicepath).name  barra contpath(k).name];
      end
   end
   indicepath = indicepath + 1;
end
listapath(:).name;

if(visualizacion)
    length(listapath)
end

archivos = [];
for directorio = 1:length(listapath)
    contpath = dir(listapath(directorio).name);
    for k = 1:length(contpath)
        if(iscellstr(extension))
            for i_ext=1:length(extension)
                ext = char(extension(i_ext));
                if( ~( contpath(k).isdir==1) ) 
		    if ( length(contpath(k).name)>length(ext) & strcmp(ext,contpath(k).name((end-length(ext)+1):end)) )
                        archivos(end+1).name = [listapath(directorio).name barra contpath(k).name];
	    	    end
                end
            end
        else
            if( ~( contpath(k).isdir==1) ) 
		if ( length(contpath(k).name)>length(extension) & strcmp(extension,contpath(k).name((end-length(extension)+1):end)) )
	            archivos(end+1).name = [listapath(directorio).name barra contpath(k).name];
		end
            end
        end
    end
end

% Se agrega al struct el nombre de referencia.
for i = 1 : length(archivos)
    path_wav = archivos(i).name;
    inds_barra = find(path_wav == barra);
    inds_point = find(path_wav == '.');
    archivos(i).clean_file_name = path_wav(inds_barra(end)+1:inds_point(end)-1);
end


end
