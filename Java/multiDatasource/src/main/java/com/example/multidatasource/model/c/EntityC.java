
package com.example.multidatasource.model.c;

import lombok.Data;
import javax.persistence.*;

@Data
@Entity
@Table(name = "existing_table")
public class EntityC {
    @Id
    private Long id;
    
    @Column(name = "existing_column1")
    private String field1;
    
    @Column(name = "existing_column2")
    private String field2;
}
