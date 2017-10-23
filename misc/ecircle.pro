; ecircle.pro,  started 7-23-98 by ERW
;
; kludged from circarray()


function ecircle,diameter,radius=radius,scale=scale,center=center, $
	dimension=dimension,truncate=truncate
; returns an array, size equal to "array" variable, with value 1
; everywhere within a circle of diameter of the array size.  Circle
; is at center of array.
;
; 'radius' sets a radius different from the default radius (half the array size)
; 'scale' lets you set a vector for the relative rescaling of x,y,z...
;    ----> [1,1,2] squashes the sphere in the z-direction
; 'center' overrides the default center of the circle
; 'dimension' lets you choose the dimension...
; /truncate to return only the parallelpiped containing data (if you
;    are using 'scale' you might want this.)

; kludgy below:
array=fltarr(diameter,diameter)
if (keyword_set(dimension)) then begin
	if (dimension eq 1) then array=fltarr(diameter)
	if (dimension eq 3) then array=fltarr(diameter,diameter,diameter)
endif

s=size(array)
result=array*0

case s(0) of

	1: begin
		sx=s(1)
		minsize = sx
		cx=(sx-1)*0.5
		if keyword_set(center) then cx=center(0)
		if keyword_set(radius) then begin
			irad=radius
		endif else begin
			irad=minsize/2
		endelse
		irad = irad*irad
		for k=0,sx-1 do begin
			rad3 = 1.0*(cx-k)*(cx-k)
			if keyword_set(scale) then rad3 = rad3*scale(0)*scale(0)
			result(k) = (rad3 le irad)
		endfor
	end

	2: begin
		sx=s(1) & sy=s(2)
		minsize = (sx < sy)
		cx=(sx-1)*0.5 & cy=(sy-1)*0.5
		if keyword_set(center) then begin
			cx=center(0)
			cy=center(1)
		endif
		if keyword_set(radius) then begin
			irad=radius
		endif else begin
			irad=minsize/2
		endelse
		irad=irad*1.0
		irad = irad*irad*1.0
		for j=0,sy-1 do begin
			rad2 = (cy-j*1.0)*(cy-j*1.0)*1.0
			if keyword_set(scale) then rad2 = 1.0*rad2*scale(1)*scale(1)
			for k=0,sx-1 do begin
				if keyword_set(scale) then begin
					rad3 = rad2 + 1.0*(cx-k*1.0)*(cx-k*1.0)*scale(0)*scale(0)
				endif else begin
					rad3 = rad2 + 1.0*(cx-k*1.0)*(cx-k*1.0)
				endelse
				result(k,j) = (rad3 le irad)
			endfor
		endfor
	end

	3: begin
		sx=s(1) & sy=s(2) & sz=s(3)
		minsize = (sx < sy < sz)
		cx=(sx-1)*0.5 & cy=(sy-1)*0.5 & cz=(sz-1)*0.5
		if keyword_set(center) then begin
			cx=center(0)
			cy=center(1)
			cz=center(2)
		endif
		if keyword_set(radius) then begin
			irad=radius
		endif else begin
			irad=minsize/2
		endelse
		irad = irad*irad
		for i=0,sz-1 do begin
			rad1=1.0*(cz-i)*(cz-i)
			if keyword_set(scale) then rad1 = rad1*scale(2)*scale(2)
			for j=0,sy-1 do begin
				rad2 = rad1 + (cy-j)*(cy-j)
				if keyword_set(scale) then begin
					rad2 = rad1 + (cy-j)*(cy-j)*scale(1)*scale(1)
				endif else begin
					rad2 = rad1 + (cy-j)*(cy-j)
				endelse
				for k=0,sx-1 do begin
					if keyword_set(scale) then begin
						rad3 = rad2 + (cx-k)*(cx-k)*scale(0)*scale(0)
					endif else begin
						rad3 = rad2 + (cx-k)*(cx-k)
					endelse
					result(k,j,i) = (rad3 le irad)
				endfor
			endfor
		endfor
	end

else: message,'dimension of array is too big, sorry!',/inf

endcase

if (keyword_set(truncate)) then begin
	s=size(result)
	w=where(result gt 0)
	zz=w / (s(1)*s(2))
	yy=(w-(zz*s(1)*s(2))) / s(1)
	xx=(w-(zz*s(1)*s(2))) - (yy*s(1))
	bx=max(xx,min=ax)
	by=max(yy,min=ay)
	bz=max(zz,min=az)
	result=result(ax:bx,ay:by,az:bz)
endif

return,result

end


